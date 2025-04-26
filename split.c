/**
 * @file split_merged.c
 * @brief Processes PGM images to detect and extract vertical tracks.
 *
 * This program analyzes a grayscale PGM (P5 format) image to find vertical segments
 * (tracks) based on brightness thresholds and saves each detected track as a
 * separate PGM file. It supports two main modes for edge detection, switchable
 * via the '-e' flag:
 *
 * 1. V1 Mode (default): Uses a static left edge detection for each segment based
 * on the first column containing a dark pixel. Column averages are smoothed
 * using the '-s' option.
 * 2. V2 Mode (-e): Uses dynamic row-by-row left edge detection with horizontal
 * smoothing ('-s'), search delta ('-d'), minimum run length ('-k'), and
 * median filtering ('-M') for edge stabilization. Initial column average
 * smoothing is controlled by '-S'.
 *
 * Common options include input file (-i), output directory (-m), brightness
 * threshold (-t), forced output width (-w), and debug (-b).
 */

#define _POSIX_C_SOURCE 200809L // For getopt, strdup, basename, snprintf
#include <stdio.h>
#include <stdlib.h> // For qsort, exit, atoi, round, strtol, malloc, realloc, free
#include <string.h> // For strncmp, memcpy, memset, strrchr, strlen, strerror, strncpy
#include <math.h>   // For round
#include <unistd.h> // For getopt
#include <errno.h>
#include <ctype.h>   // For isprint, tolower, isspace
#include <sys/stat.h> // For mkdir, stat
#include <libgen.h>   // For basename
#include <limits.h> // For INT_MAX

// For mkdir in Windows
#ifdef _WIN32
#include <direct.h>
#define MKDIR(path, mode) _mkdir(path)
#else
#define MKDIR(path, mode) mkdir(path, mode)
#endif

// --- Constants ---
#define MAX_PATH_LEN 512
#define MAX_LINE_LEN 256
#define PGM_MAXVAL_SUPPORTED 255
#define WIDTH_REDUCTION_FACTOR 0.95 // Factor for auto-width calculation

// Default parameter values
#define DEFAULT_V1_SMOOTH_WINDOW 1          // Default for -s in V1 mode
#define DEFAULT_V2_INITIAL_SMOOTH_WINDOW 1  // Default for -S in V2 mode
#define DEFAULT_V2_ROW_SMOOTH_WINDOW 1      // Default for -s in V2 mode
#define DEFAULT_PROCESSING_WIDTH 0          // Default for -w (auto)
#define DEFAULT_V2_SEARCH_DELTA_DIVIDER 2   // Default calculation for -d
#define DEFAULT_V2_MIN_RUN_LENGTH 1         // Default for -k
#define DEFAULT_V2_MEDIAN_WINDOW 3          // Default for -M

// --- Structures ---

/**
 * @brief Configuration parameters from command line options.
 */
typedef struct {
    char input_filename[MAX_PATH_LEN];      // -i: Input PGM file path
    char output_m_directory[MAX_PATH_LEN];  // -m: Output directory path
    int threshold_abs;                      // -t: Absolute brightness threshold
    int initial_smooth_window;              // Effective initial smoothing window (-S or -s)
    int row_smooth_window;                  // -s: Row smoothing window (V2 only)
    int median_window;                      // -M: Median filter window (V2 only)
    int processing_width;                   // -w: Forced output width (0=auto)
    int search_delta;                       // -d: Edge search delta (V2 only)
    int min_run_length;                     // -k: Min dark run length (V2 only)
    int debug;                              // -b: Debug flag
    int use_v2_logic;                       // -e: Flag to enable V2 logic

    // Calculated values used during processing
    int calculated_avg_width;               // Average width between track starts
    int final_processing_width;             // Final width used for output images
} Config;

/**
 * @brief Structure to hold PGM image data.
 */
typedef struct {
    unsigned char *data; // Pixel data (grayscale)
    int width;           // Image width
    int height;          // Image height
    int maxval;          // Maximum pixel value (from PGM header)
} Image;

// --- Function Prototypes ---

// Helpers
void exit_error(const char *message, const char *details);
void* safe_malloc(size_t size);
void* safe_realloc(void* ptr, size_t size);
int int_compare(const void *a, const void *b);
long safe_strtol(const char *str, const char *arg_name);
void get_pgm_basename(const char* filepath, char* out_basename, size_t max_len);
void print_usage(const char *prog_name);

// PGM I/O
Image* read_pgm_p5(const char *filename);
int write_pgm_p5(const char *filename, const unsigned char *data, int width, int height, int maxval);

// Image Analysis
int is_blank_row(const unsigned char *data, int img_width, int row_index, int total_height, int col_start, int col_end, int threshold);
double calculate_column_average_range(const unsigned char *data, int img_width, int col, int row_start, int row_end);
void smooth_1d_array(const double *data_in, double *data_out, int n, int window_size);
int column_has_dark_pixel(const unsigned char* track_data, int track_w, int track_h, int col, int threshold); // V1 specific helper
void smooth_horizontal_row(const unsigned char *img_data, int img_width, int img_height, int row, int col_start, int col_end, int window_size, double *out_smoothed_row); // V2 specific helper
int find_left_edge_in_row(const unsigned char *img_data, int img_width, int img_height, int r_img, int predicted_col_orig, const Config *config, double *smoothed_row_buffer); // V2 specific helper

// Segment Processing
void process_segment_v1(const unsigned char* cropped_data, int cropped_w, int cropped_h, int segment_start_rel, int segment_end_rel, int original_maxval, const Config* config, int track_num, int average_width_hint, const char* input_basename);
void process_segment_v2(const Image* full_img, int global_top, int global_bottom, int segment_start_rel, int left, int track_index, const Config* config, const char* input_basename);


// --- Helper Function Definitions ---

void exit_error(const char *message, const char *details) {
    fprintf(stderr, "ERROR: %s%s%s\n",
            message ? message : "Unknown error",
            details ? ": " : "",
            details ? details : "");
    exit(EXIT_FAILURE);
}

void* safe_malloc(size_t size) {
    if (size == 0) size = 1;
    void *ptr = malloc(size);
    if (!ptr) {
        exit_error("Memory allocation failed", strerror(errno));
    }
    return ptr;
}

void* safe_realloc(void* ptr, size_t size) {
    if (size == 0) {
        free(ptr);
        return NULL;
    }
    void *new_ptr = realloc(ptr, size);
    if (!new_ptr) {
        exit_error("Memory reallocation failed", strerror(errno));
    }
    return new_ptr;
}

int int_compare(const void *a, const void *b) {
    int int_a = *((const int*)a);
    int int_b = *((const int*)b);
    if (int_a < int_b) return -1;
    if (int_a > int_b) return 1;
    return 0;
}

long safe_strtol(const char *str, const char *arg_name) {
    char *endptr;
    errno = 0;
    long val = strtol(str, &endptr, 10);

    if (errno == ERANGE) {
        char err_details[100];
        snprintf(err_details, sizeof(err_details), "Value '%s' for %s out of range.", str, arg_name);
        exit_error("Invalid number format", err_details);
    }
    if (errno != 0 && val == 0) {
         char err_details[100];
         snprintf(err_details, sizeof(err_details), "Invalid value '%s' for %s.", str, arg_name);
         exit_error("Invalid number format", err_details);
    }
    if (endptr == str) {
        char err_details[100];
        snprintf(err_details, sizeof(err_details), "No digits found in '%s' for %s.", str, arg_name);
        exit_error("Invalid number format", err_details);
    }
    while (isspace((unsigned char)*endptr)) {
        endptr++;
    }
    if (*endptr != '\0') {
        char err_details[100];
        snprintf(err_details, sizeof(err_details), "Trailing characters '%s' after number in '%s' for %s.", endptr, str, arg_name);
        exit_error("Invalid number format", err_details);
    }
    return val;
}

void get_pgm_basename(const char* filepath, char* out_basename, size_t max_len) {
    if (!filepath || !out_basename || max_len == 0) return;
    char filepath_copy[MAX_PATH_LEN];
    strncpy(filepath_copy, filepath, sizeof(filepath_copy) - 1);
    filepath_copy[sizeof(filepath_copy) - 1] = '\0';

    char *base = basename(filepath_copy);
    strncpy(out_basename, base, max_len - 1);
    out_basename[max_len - 1] = '\0';

    char *dot = strrchr(out_basename, '.');
    if (dot != NULL && dot != out_basename && strlen(dot) == 4 &&
        (tolower((unsigned char)dot[1]) == 'p') &&
        (tolower((unsigned char)dot[2]) == 'g') &&
        (tolower((unsigned char)dot[3]) == 'm'))
    {
        *dot = '\0';
    }
}

void print_usage(const char *prog_name) {
     fprintf(stderr, "\nUsage: %s -i <in.pgm> -m <out_dir> -t <thresh> [-e] [options]\n", prog_name);
     fprintf(stderr, "Required:\n");
     fprintf(stderr, "  -i <input.pgm>    : Input PGM file (P5 format, maxval <= %d)\n", PGM_MAXVAL_SUPPORTED);
     fprintf(stderr, "  -m <pgm_dir>      : Output directory for segments\n");
     fprintf(stderr, "  -t <threshold>    : Absolute brightness threshold [0-%d]\n", PGM_MAXVAL_SUPPORTED);
     fprintf(stderr, "Mode Selection:\n");
     fprintf(stderr, "  -e                : Use V2 dynamic edge tracking logic (default: V1 static edge)\n");
     fprintf(stderr, "Common Optional:\n");
     fprintf(stderr, "  -w <width>        : Segment processing width (>=0, 0=auto ~avg*%.2f, default: %d)\n", WIDTH_REDUCTION_FACTOR, DEFAULT_PROCESSING_WIDTH);
     fprintf(stderr, "  -b                : Enable debug messages\n");
     fprintf(stderr, "V1 Logic Options (used if -e is NOT set):\n");
     fprintf(stderr, "  -s <smooth_V1>    : Initial column average smoothing window (>=1, default: %d)\n", DEFAULT_V1_SMOOTH_WINDOW);
     fprintf(stderr, "V2 Logic Options (used if -e IS set):\n");
     fprintf(stderr, "  -S <smooth_I>     : Initial column average smoothing window (>=1, default: %d)\n", DEFAULT_V2_INITIAL_SMOOTH_WINDOW);
     fprintf(stderr, "  -s <smooth_R>     : Row horizontal smoothing window for edge tracking (>=1, default: %d)\n", DEFAULT_V2_ROW_SMOOTH_WINDOW);
     fprintf(stderr, "  -M <median_W>     : Median filter window size for edge tracking (odd >=3, default: %d)\n", DEFAULT_V2_MEDIAN_WINDOW);
     fprintf(stderr, "  -d <delta>        : Edge search window delta (+/-, optional, default: final_width/%d)\n", DEFAULT_V2_SEARCH_DELTA_DIVIDER);
     fprintf(stderr, "  -k <runlen>       : Min dark pixels run length to confirm edge (>=1, default: %d)\n", DEFAULT_V2_MIN_RUN_LENGTH);
}

// --- PGM I/O Function Definitions ---

Image* read_pgm_p5(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        exit_error("Cannot open input PGM file", filename);
    }

    Image *img = (Image*)safe_malloc(sizeof(Image));
    img->data = NULL;
    img->width = 0;
    img->height = 0;
    img->maxval = 0;
    char line[MAX_LINE_LEN];
    int fields_read;
    int read_maxval = -1;

    // 1. Check Magic Number "P5"
    if (!fgets(line, sizeof(line), fp) || strncmp(line, "P5", 2) != 0) {
        fclose(fp); free(img);
        exit_error("Input file is not a P5 PGM", filename);
    }

    // 2. Read Dimensions (width height), skipping comments
    fields_read = 0;
    do {
        if (!fgets(line, sizeof(line), fp)) {
            fclose(fp); free(img);
            exit_error("Could not read dimensions/maxval from PGM header", filename);
        }
        if (line[0] != '#') {
            fields_read = sscanf(line, "%d %d", &img->width, &img->height);
        }
    } while (line[0] == '#' || fields_read < 2);

    if (fields_read != 2 || img->width <= 0 || img->height <= 0) {
        fclose(fp); free(img);
        exit_error("Invalid image dimensions in PGM file", filename);
    }

    // 3. Read Maxval, skipping comments
    fields_read = 0;
    do {
        if (!fgets(line, sizeof(line), fp)) {
            fclose(fp); free(img);
            exit_error("Could not read maxval from PGM header", filename);
        }
        if (line[0] != '#') {
            fields_read = sscanf(line, "%d", &read_maxval);
        }
    } while (line[0] == '#' || fields_read < 1);

    if (fields_read != 1 || read_maxval <= 0 || read_maxval > PGM_MAXVAL_SUPPORTED) {
        fclose(fp); free(img);
        char err_details[100];
        snprintf(err_details, sizeof(err_details), "Value %d. Must be > 0 and <= %d.", read_maxval, PGM_MAXVAL_SUPPORTED);
        exit_error("Invalid or unsupported maxval in PGM file", err_details);
    }
    img->maxval = read_maxval;

    // 4. Allocate memory for pixel data
    size_t size = (size_t)img->width * img->height;
    if (size == 0) { fclose(fp); free(img); exit_error("Image size is zero based on dimensions", filename); }
    img->data = (unsigned char*)safe_malloc(size);

    // 5. Skip single whitespace character after maxval
    int c = fgetc(fp);
    if (c == EOF) {
        fclose(fp); free(img->data); free(img);
        exit_error("Unexpected EOF before pixel data", filename);
    }
    if (!isspace(c)) {
        ungetc(c, fp);
        fprintf(stderr, "Warning: Expected single whitespace after maxval in %s, found '%c'.\n", filename, c);
    }

    // 6. Read pixel data
    size_t bytes_read = fread(img->data, 1, size, fp);
    if (bytes_read != size) {
        // *** Исправленный блок обработки ошибки fread ***
        int read_error = ferror(fp); // Проверяем ошибку ПОТОКА
        int end_of_file = feof(fp);   // Проверяем конец файла
        char error_details[100] = {0}; // Буфер для strerror

        if (read_error) {
            // Копируем сообщение об ошибке ДО fclose, т.к. fclose может изменить errno
            snprintf(error_details, sizeof(error_details), "%s", strerror(errno));
        }

        // Закрываем файл и освобождаем память ПЕРЕД выходом
        fclose(fp);
        free(img->data);
        free(img);

        // Выходим с соответствующей ошибкой
        if (read_error) {
            exit_error("Failed to read pixel data", error_details);
        } else if (end_of_file) {
            exit_error("Insufficient pixel data in file (reached EOF prematurely)", filename);
        } else {
            // Эта ветка маловероятна, если fread вернул < size без ошибки и без EOF
            exit_error("Failed to read expected amount of pixel data for unknown reason", filename);
        }
        // *** Конец исправленного блока ***
    }

    // Если чтение прошло успешно, закрываем файл
    if (fclose(fp) != 0) {
         // Ошибка закрытия файла - не критично, но сообщим
         fprintf(stderr, "Warning: Failed to close input PGM file %s properly (%s).\n", filename, strerror(errno));
    }

    return img; // Возвращаем успешно прочитанное изображение
}


int write_pgm_p5(const char *filename, const unsigned char *data, int width, int height, int maxval) {
    if (width <= 0 || height <= 0 || !data) {
        fprintf(stderr, "Warning: Invalid image data for writing PGM (%dx%d) to %s\n", width, height, filename);
        return -1;
    }
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Warning: Cannot create output PGM file %s (%s).\n", filename, strerror(errno));
        return -1;
    }

    int write_maxval = (maxval > PGM_MAXVAL_SUPPORTED) ? PGM_MAXVAL_SUPPORTED : maxval;
    fprintf(fp, "P5\n%d %d\n%d\n", width, height, write_maxval);

    size_t size = (size_t)width * height;
    if (fwrite(data, 1, size, fp) != size) {
        fprintf(stderr, "Warning: Failed to write all pixel data to %s (%s).\n", filename, strerror(errno));
        fclose(fp);
        return -1;
    }

    if (fclose(fp) != 0) {
         fprintf(stderr, "Warning: Failed to close output PGM file %s properly (%s).\n", filename, strerror(errno));
         return -1;
    }
    return 0;
}


// --- Image Analysis Function Definitions ---

int is_blank_row(const unsigned char *data, int img_width, int row_index, int total_height, int col_start, int col_end, int threshold) {
     if (row_index < 0 || row_index >= total_height) return 1;
     col_start = (col_start < 0) ? 0 : col_start;
     col_end = (col_end >= img_width) ? img_width - 1 : col_end;
     if (col_start > col_end) return 1;

     size_t row_offset = (size_t)row_index * img_width;
     for (int c = col_start; c <= col_end; c++) {
         if (data[row_offset + c] < threshold) {
             return 0;
         }
     }
     return 1;
}

double calculate_column_average_range(const unsigned char *data, int img_width, int col, int row_start, int row_end) {
    if (col < 0 || col >= img_width || row_start < 0 || row_start > row_end) {
         return (double)PGM_MAXVAL_SUPPORTED + 1.0;
    }
    int height = row_end - row_start + 1;
    if (height <= 0) {
         return (double)PGM_MAXVAL_SUPPORTED + 1.0;
    }

    long long sum = 0;
    for (int r = row_start; r <= row_end; r++) {
        sum += data[(size_t)r * img_width + col];
    }
    return (double)sum / height;
}

void smooth_1d_array(const double *data_in, double *data_out, int n, int window_size) {
    if (n <= 0) return;
    if (window_size <= 1) {
        memcpy(data_out, data_in, n * sizeof(double));
        return;
    }

    int half_window_floor = window_size / 2;
    int half_window_ceil = (window_size - 1) / 2;

    for (int i = 0; i < n; ++i) {
        double sum = 0;
        int count = 0;
        int k_start = (i - half_window_floor < 0) ? 0 : i - half_window_floor;
        int k_end = (i + half_window_ceil >= n) ? n - 1 : i + half_window_ceil;

        for (int k = k_start; k <= k_end; ++k) {
            sum += data_in[k];
            count++;
        }
        data_out[i] = (count > 0) ? sum / count : data_in[i];
    }
}

int column_has_dark_pixel(const unsigned char* track_data, int track_w, int track_h, int col, int threshold) {
    if (col < 0 || col >= track_w || track_h <= 0) {
        return 0;
    }
    for (int r = 0; r < track_h; r++) {
        if (track_data[(size_t)r * track_w + col] < threshold) {
            return 1;
        }
    }
    return 0;
}

void smooth_horizontal_row(const unsigned char *img_data, int img_width, int img_height, int row, int col_start, int col_end, int window_size, double *out_smoothed_row) {
     if (row < 0 || row >= img_height) return;
     col_start = (col_start < 0) ? 0 : col_start;
     col_end = (col_end >= img_width) ? img_width - 1 : col_end;
     if (col_start > col_end) return;

     int out_len = col_end - col_start + 1;

     if (window_size <= 1) {
         size_t row_offset = (size_t)row * img_width;
         for(int i=0; i<out_len; ++i) {
             out_smoothed_row[i] = (double)img_data[row_offset + col_start + i];
         }
         return;
     }

     int half_window_floor = window_size / 2;
     int half_window_ceil = (window_size - 1) / 2;
     size_t row_offset = (size_t)row * img_width;

     for (int c_idx = 0; c_idx < out_len; ++c_idx) {
         int c_orig = col_start + c_idx;
         double sum = 0;
         int count = 0;
         int k_start_orig = c_orig - half_window_floor;
         int k_end_orig = c_orig + half_window_ceil;

         for (int k_orig = k_start_orig; k_orig <= k_end_orig; ++k_orig) {
             if (k_orig >= 0 && k_orig < img_width) {
                 sum += img_data[row_offset + k_orig];
                 count++;
             }
         }
         out_smoothed_row[c_idx] = (count > 0) ? sum / count : (double)img_data[row_offset + c_orig];
     }
}

int find_left_edge_in_row(
    const unsigned char *img_data, int img_width, int img_height,
    int r_img, int predicted_col_orig,
    const Config *config, double *smoothed_row_buffer
) {
    if (r_img < 0 || r_img >= img_height) return -1;

    int search_start_orig = predicted_col_orig - config->search_delta;
    int search_end_orig = predicted_col_orig + config->search_delta;
    search_start_orig = (search_start_orig < 0) ? 0 : search_start_orig;
    search_end_orig = (search_end_orig >= img_width) ? img_width - 1 : search_end_orig; // Исправлено здесь
    if (search_start_orig > search_end_orig) return -1;

    int search_len = search_end_orig - search_start_orig + 1;

    smooth_horizontal_row(img_data, img_width, img_height, r_img, search_start_orig, search_end_orig, config->row_smooth_window, smoothed_row_buffer);

    int found_edge_orig = -1;
    double prev_smoothed_val;

    if (search_start_orig > 0) {
        double temp_buf[1];
        smooth_horizontal_row(img_data, img_width, img_height, r_img, search_start_orig - 1, search_start_orig - 1, config->row_smooth_window, temp_buf);
        prev_smoothed_val = temp_buf[0];
    } else {
        prev_smoothed_val = (double)PGM_MAXVAL_SUPPORTED + 1.0;
    }

    for (int i = 0; i < search_len; ++i) {
        double current_smoothed_val = smoothed_row_buffer[i];
        int current_col_orig = search_start_orig + i;

        if (current_smoothed_val < config->threshold_abs && prev_smoothed_val >= config->threshold_abs) {
            if (config->min_run_length > 1) {
                int dark_run_count = 0;
                size_t row_offset = (size_t)r_img * img_width;
                for (int k = 0; k < config->min_run_length; ++k) {
                    int check_col_orig = current_col_orig + k;
                    if (check_col_orig >= img_width) break;
                    if (img_data[row_offset + check_col_orig] < config->threshold_abs) {
                        dark_run_count++;
                    } else {
                        break;
                    }
                }
                if (dark_run_count >= config->min_run_length) {
                    found_edge_orig = current_col_orig;
                    break;
                }
            } else {
                found_edge_orig = current_col_orig;
                break;
            }
        }
        prev_smoothed_val = current_smoothed_val;
    }
    return found_edge_orig;
}


// --- Segment Processing Function Definitions ---

void process_segment_v1(
    const unsigned char* cropped_data, int cropped_w, int cropped_h,
    int segment_start_rel, int segment_end_rel,
    int original_maxval, const Config* config, int track_num,
    int average_width_hint, const char* input_basename
) {
    if (config->debug) {
        printf(" V1: Processing Track %d (Cols %d-%d relative, Height %d)...\n",
               track_num, segment_start_rel, segment_end_rel, cropped_h);
    }

    int initial_segment_w = segment_end_rel - segment_start_rel + 1;
    if (initial_segment_w <= 0 || cropped_h <= 0) {
        fprintf(stderr,"Warning (V1): Invalid initial segment dimensions for track %d (%dx%d). Skipping.\n",
                track_num, initial_segment_w, cropped_h);
        return;
    }

    unsigned char* segment_data = (unsigned char*)safe_malloc((size_t)initial_segment_w * cropped_h);
    for (int r = 0; r < cropped_h; r++) {
        memcpy(segment_data + (size_t)r * initial_segment_w,
               cropped_data + (size_t)r * cropped_w + segment_start_rel,
               initial_segment_w);
    }

    int left_edge_rel = 0;
    while(left_edge_rel < initial_segment_w &&
          !column_has_dark_pixel(segment_data, initial_segment_w, cropped_h, left_edge_rel, config->threshold_abs)) {
        left_edge_rel++;
    }

    int right_boundary_rel = initial_segment_w - 1;
    if (left_edge_rel > right_boundary_rel) {
         if (config->debug) printf(" Skipped track %d (V1): segment empty after left-crop check (col %d relative).\n", track_num, left_edge_rel);
         free(segment_data); return;
    }
    if (config->debug > 1) printf("  Track %d (V1): Refined left edge (relative to segment start): %d\n", track_num, left_edge_rel);

    int top_crop = 0;
    int bottom_crop = cropped_h - 1;
    while (top_crop <= bottom_crop && is_blank_row(segment_data, initial_segment_w, top_crop, cropped_h, left_edge_rel, right_boundary_rel, config->threshold_abs)) {
        top_crop++;
    }
    while (bottom_crop >= top_crop && is_blank_row(segment_data, initial_segment_w, bottom_crop, cropped_h, left_edge_rel, right_boundary_rel, config->threshold_abs)) {
        bottom_crop--;
    }

    int final_h = bottom_crop - top_crop + 1;
    if (final_h <= 0) {
        if (config->debug) printf(" Skipped track %d (V1): segment empty after vertical crop (rows %d-%d relative).\n", track_num, top_crop, bottom_crop);
        free(segment_data); return;
    }
    if (config->debug > 1) printf("  Track %d (V1): Vertically cropped (relative rows): %d-%d (height %d)\n", track_num, top_crop, bottom_crop, final_h);

    int final_w;
    char width_source[50];
    int actual_content_width_rel = right_boundary_rel - left_edge_rel + 1;

    if (config->processing_width > 0) {
        final_w = config->processing_width;
        snprintf(width_source, sizeof(width_source), "forced by -w");
    } else if (average_width_hint > 0) {
        int reduced_width = (int)round((double)average_width_hint * WIDTH_REDUCTION_FACTOR);
        final_w = (reduced_width <= 0) ? 1 : reduced_width;
        snprintf(width_source, sizeof(width_source), "avg*%.2f", WIDTH_REDUCTION_FACTOR);
    } else {
        final_w = actual_content_width_rel;
        snprintf(width_source, sizeof(width_source), "content width");
    }
    if (final_w <= 0) final_w = 1;

    unsigned char *final_track_data = (unsigned char*)safe_malloc((size_t)final_w * final_h);
    int bg_color = (original_maxval > PGM_MAXVAL_SUPPORTED) ? PGM_MAXVAL_SUPPORTED : original_maxval;
    memset(final_track_data, bg_color, (size_t)final_w * final_h);

    int copy_w = (final_w < actual_content_width_rel) ? final_w : actual_content_width_rel;
    if (copy_w > 0) {
        for (int r = 0; r < final_h; r++) {
            unsigned char* source_ptr = segment_data + (size_t)(top_crop + r) * initial_segment_w + left_edge_rel;
            unsigned char* dest_ptr = final_track_data + (size_t)r * final_w;
            memcpy(dest_ptr, source_ptr, copy_w);
        }
    }

    char fname[MAX_PATH_LEN * 2];
    int chars_written = snprintf(fname, sizeof(fname), "%s/%s_%04d.pgm",
                                 config->output_m_directory, input_basename, track_num);

    if (chars_written < 0 || chars_written >= (int)sizeof(fname)) {
         fprintf(stderr, "Warning (V1): Failed to format filename for track %d. Skipping save.\n", track_num);
    } else {
        if (write_pgm_p5(fname, final_track_data, final_w, final_h, original_maxval) == 0) {
             printf("Saved track (V1 logic): %s (width %d - %s)\n", fname, final_w, width_source);
        }
    }

    free(segment_data);
    free(final_track_data);
}


void process_segment_v2(
    const Image* full_img,
    int global_top, int global_bottom,
    int segment_start_rel,
    int left,
    int track_index, const Config* config,
    const char* input_basename
) {
     int predicted_start_col_orig = left + segment_start_rel;
     if (config->debug) {
         printf(" V2: Processing Track %d (predicted start col: %d absolute)...\n",
                track_index, predicted_start_col_orig);
     }

    int vc_col_start = predicted_start_col_orig - config->search_delta;
    int vc_col_end   = predicted_start_col_orig + config->final_processing_width + config->search_delta;
    vc_col_start = (vc_col_start < 0) ? 0 : vc_col_start;
    vc_col_end = (vc_col_end >= full_img->width) ? full_img->width - 1 : vc_col_end;

    int segment_top_local = global_top;
    int segment_bottom_local = global_bottom;
    while (segment_top_local <= segment_bottom_local && is_blank_row(full_img->data, full_img->width, segment_top_local, full_img->height, vc_col_start, vc_col_end, config->threshold_abs)) segment_top_local++;
    while (segment_bottom_local >= segment_top_local && is_blank_row(full_img->data, full_img->width, segment_bottom_local, full_img->height, vc_col_start, vc_col_end, config->threshold_abs)) segment_bottom_local--;

    int segment_h_local = segment_bottom_local - segment_top_local + 1;
    if (segment_h_local <= 0) {
        fprintf(stderr, "Warning (V2): Track %d skipped. Height <= 0 after vertical crop (%d-%d).\n",
                track_index, segment_top_local, segment_bottom_local);
        return;
    }
    if (config->debug) {
        printf("  Track %d (V2): Vertically cropped to rows %d-%d (height %d)\n",
               track_index, segment_top_local, segment_bottom_local, segment_h_local);
    }

    unsigned char *pgm_segment_data = NULL;
    char pgm_filename[MAX_PATH_LEN * 2];
    pgm_segment_data = (unsigned char*)safe_malloc((size_t)config->final_processing_width * segment_h_local);
    int bg_color = (full_img->maxval > PGM_MAXVAL_SUPPORTED) ? PGM_MAXVAL_SUPPORTED : full_img->maxval;
    memset(pgm_segment_data, bg_color, (size_t)config->final_processing_width * segment_h_local);

    int snprintf_ret = snprintf(pgm_filename, sizeof(pgm_filename), "%s/%s_%04d.pgm",
                                config->output_m_directory, input_basename, track_index);
    if (snprintf_ret < 0 || snprintf_ret >= (int)sizeof(pgm_filename)) {
        fprintf(stderr, "Warning (V2): Failed to create filename for PGM track %d. PGM saving disabled.\n", track_index);
        free(pgm_segment_data);
        pgm_segment_data = NULL;
    }

    int *edge_history = (int*)safe_malloc(config->median_window * sizeof(int));
    int *sort_buffer = (int*)safe_malloc(config->median_window * sizeof(int));
    int history_count = 0;
    int history_idx = 0;
    int current_median_edge = predicted_start_col_orig;

    int max_search_len = (config->search_delta * 2) + 1;
    double *smoothed_row_buffer = (double*)safe_malloc(max_search_len * sizeof(double));

    for (int r_img = segment_top_local; r_img <= segment_bottom_local; ++r_img) {
        int r_rel = r_img - segment_top_local;

        int raw_edge = find_left_edge_in_row(full_img->data, full_img->width, full_img->height,
                                           r_img, current_median_edge, config, smoothed_row_buffer);

        int edge_to_store = (raw_edge != -1) ? raw_edge : current_median_edge;
        edge_history[history_idx] = edge_to_store;
        history_idx = (history_idx + 1) % config->median_window;
        if (history_count < config->median_window) history_count++;

        if (history_count > 0) {
            int sort_idx = 0;
            int start_copy_idx = (history_count == config->median_window) ? history_idx : 0;
            for(int k=0; k<history_count; ++k) {
                int current_copy_idx = (start_copy_idx + k) % config->median_window;
                sort_buffer[sort_idx++] = edge_history[current_copy_idx];
            }
            qsort(sort_buffer, history_count, sizeof(int), int_compare);
            current_median_edge = sort_buffer[history_count / 2];
        }

        if (config->debug > 1) {
            printf("  Track %d, Row %d: Raw=%d, Stored=%d, Median(%d pts)=%d\n",
                   track_index, r_img, raw_edge, edge_to_store, history_count, current_median_edge);
        }

        if (pgm_segment_data) {
            int win_start_col_orig = current_median_edge;
            unsigned char *dest_row_ptr = pgm_segment_data + (size_t)r_rel * config->final_processing_width;
            size_t src_row_offset = (size_t)r_img * full_img->width;

            for (int out_col = 0; out_col < config->final_processing_width; ++out_col) {
                int src_col = win_start_col_orig + out_col;
                if (src_col >= 0 && src_col < full_img->width) {
                    dest_row_ptr[out_col] = full_img->data[src_row_offset + src_col];
                } else {
                    dest_row_ptr[out_col] = (unsigned char)bg_color;
                }
            }
        }
    }

    if (pgm_segment_data) {
        if (write_pgm_p5(pgm_filename, pgm_segment_data, config->final_processing_width, segment_h_local, full_img->maxval) == 0) {
             printf("Saved track (V2 logic): %s (%dx%d)\n", pgm_filename, config->final_processing_width, segment_h_local);
        }
        free(pgm_segment_data);
    }
    free(edge_history);
    free(sort_buffer);
    free(smoothed_row_buffer);
}


// --- Main Function ---
int main(int argc, char *argv[]) {
    Config config;
    long temp_long;

    // Set defaults
    config.input_filename[0] = '\0';
    config.output_m_directory[0] = '\0';
    config.threshold_abs = -1;
    config.initial_smooth_window = DEFAULT_V2_INITIAL_SMOOTH_WINDOW;
    config.row_smooth_window = DEFAULT_V2_ROW_SMOOTH_WINDOW;
    config.median_window = DEFAULT_V2_MEDIAN_WINDOW;
    config.processing_width = DEFAULT_PROCESSING_WIDTH;
    config.search_delta = -1;
    config.min_run_length = DEFAULT_V2_MIN_RUN_LENGTH;
    config.debug = 0;
    config.use_v2_logic = 0;
    config.calculated_avg_width = 0;
    config.final_processing_width = 0;

    // Parse arguments
    int opt;
    // Check if no arguments were provided
    if (argc == 1) {
        print_usage(argv[0]); // Показать справку
        return EXIT_SUCCESS;  // Завершить программу успешно
    }
    while ((opt = getopt(argc, argv, "i:m:t:S:s:M:w:d:k:be")) != -1) {
        switch (opt) {
            case 'i': strncpy(config.input_filename, optarg, MAX_PATH_LEN - 1); config.input_filename[MAX_PATH_LEN - 1] = '\0'; break;
            case 'm': strncpy(config.output_m_directory, optarg, MAX_PATH_LEN - 1); config.output_m_directory[MAX_PATH_LEN - 1] = '\0'; break;
            case 't':
                temp_long = safe_strtol(optarg, "-t");
                if (temp_long < 0 || temp_long > PGM_MAXVAL_SUPPORTED) exit_error("Threshold (-t) must be between 0 and PGM_MAXVAL_SUPPORTED", optarg);
                config.threshold_abs = (int)temp_long;
                break;
            case 'S':
                temp_long = safe_strtol(optarg, "-S");
                if (temp_long < 1 || temp_long > INT_MAX) exit_error("Initial smooth window (-S) must be >= 1", optarg);
                config.initial_smooth_window = (int)temp_long;
                break;
            case 's':
                temp_long = safe_strtol(optarg, "-s");
                if (temp_long < 1 || temp_long > INT_MAX) exit_error("Smooth window (-s) must be >= 1", optarg);
                config.row_smooth_window = (int)temp_long;
                break;
            case 'M':
                temp_long = safe_strtol(optarg, "-M");
                if (temp_long < 3 || temp_long > INT_MAX || temp_long % 2 == 0) exit_error("Median filter window (-M) must be an odd integer >= 3", optarg);
                config.median_window = (int)temp_long;
                break;
            case 'w':
                temp_long = safe_strtol(optarg, "-w");
                if (temp_long < 0 || temp_long > INT_MAX) exit_error("Processing width (-w) cannot be negative", optarg);
                config.processing_width = (int)temp_long;
                break;
            case 'd':
                temp_long = safe_strtol(optarg, "-d");
                if (temp_long <= 0 || temp_long > INT_MAX / 2) exit_error("Search delta (-d) must be positive and reasonable", optarg);
                config.search_delta = (int)temp_long;
                break;
            case 'k':
                temp_long = safe_strtol(optarg, "-k");
                if (temp_long < 1 || temp_long > INT_MAX) exit_error("Min run length (-k) must be >= 1", optarg);
                config.min_run_length = (int)temp_long;
                break;
            case 'b': config.debug = 1; break;
            case 'e': config.use_v2_logic = 1; break;
            case '?':
                 print_usage(argv[0]);
                 return EXIT_FAILURE;
            default:
                 abort();
        }
    }

    // Validate required args
    if (strlen(config.input_filename) == 0) exit_error("Input PGM filename (-i) is required", NULL);
    if (strlen(config.output_m_directory) == 0) exit_error("Output PGM directory (-m) is required", NULL);
    if (config.threshold_abs == -1) exit_error("Parameter -t <threshold> is required", NULL);
    if (optind < argc) {
        fprintf(stderr, "Error: Unexpected non-option arguments found: ");
        while (optind < argc) fprintf(stderr, "%s ", argv[optind++]);
        fprintf(stderr, "\n");
        exit_error("Invalid arguments", NULL);
    }

    // Adjust smoothing parameter based on mode
    if (!config.use_v2_logic) { // V1 Mode
        config.initial_smooth_window = config.row_smooth_window;
        if (config.initial_smooth_window == DEFAULT_V2_ROW_SMOOTH_WINDOW) {
             config.initial_smooth_window = DEFAULT_V1_SMOOTH_WINDOW;
        }
        printf("Running in V1 mode (static edge detection).\n");
        if (config.median_window != DEFAULT_V2_MEDIAN_WINDOW) fprintf(stderr, "Warning: Option -M ignored in V1 mode.\n");
        if (config.search_delta > 0) fprintf(stderr, "Warning: Option -d ignored in V1 mode.\n");
        if (config.min_run_length != DEFAULT_V2_MIN_RUN_LENGTH) fprintf(stderr, "Warning: Option -k ignored in V1 mode.\n");
        // Check if -S was provided but ignored
        if (config.initial_smooth_window != DEFAULT_V2_INITIAL_SMOOTH_WINDOW && config.row_smooth_window == DEFAULT_V2_ROW_SMOOTH_WINDOW) {
             // This condition is tricky, means -S was set but -s was not, and we are in V1 mode.
             // The current logic uses the default V1 value in this case. Let's keep it simple.
        }
    } else { // V2 Mode
        printf("Running in V2 mode (dynamic edge tracking with median filter).\n");
    }

    // Get input basename
    char input_basename[MAX_PATH_LEN];
    get_pgm_basename(config.input_filename, input_basename, sizeof(input_basename));

    // Print effective parameters
    printf("--- Effective Parameters ---\n");
    printf("Input PGM (-i)      : %s (basename: %s)\n", config.input_filename, input_basename);
    printf("Output PGM Dir (-m) : %s\n", config.output_m_directory);
    printf("Threshold (-t)      : %d\n", config.threshold_abs);
    printf("Processing Width (-w): %d (0=auto ~avg*%.2f)\n", config.processing_width, WIDTH_REDUCTION_FACTOR);
    printf("Initial Smooth      : %d (from %s)\n", config.initial_smooth_window, config.use_v2_logic ? "-S" : "-s");
    if(config.use_v2_logic) {
        printf("Row Smooth (-s)     : %d\n", config.row_smooth_window);
        printf("Median Filter (-M)  : %d\n", config.median_window);
        printf("Min Run Length (-k) : %d\n", config.min_run_length);
    }
    printf("Debug (-b)          : %s\n", config.debug ? "Enabled" : "Disabled");
    printf("--------------------------\n");

    // Create output directory
    struct stat st = {0};
    if (stat(config.output_m_directory, &st) == -1) {
        printf("Creating directory: %s\n", config.output_m_directory);
        int mkdir_ret = MKDIR(config.output_m_directory, 0755);
        if (mkdir_ret != 0) exit_error("Could not create directory specified by -m", strerror(errno));
    } else if (!S_ISDIR(st.st_mode)) exit_error("Path specified by -m exists but is not a directory", config.output_m_directory);

    // ================================================================
    // Phase 1: Initial Analysis
    // ================================================================
    printf("\nPhase 1: Initial Analysis...\n");
    printf("Reading input PGM: %s...\n", config.input_filename);
    Image *img = read_pgm_p5(config.input_filename);
    printf("Read image: %d x %d, maxval=%d\n", img->width, img->height, img->maxval);

    printf("Finding overall content boundaries...\n");
    double *col_averages_full = (double*)safe_malloc(img->width * sizeof(double));
    for (int c = 0; c < img->width; c++) col_averages_full[c] = calculate_column_average_range(img->data, img->width, c, 0, img->height - 1);
    int left = 0; while (left < img->width && col_averages_full[left] >= config.threshold_abs) left++;
    int right = img->width - 1; while (right >= left && col_averages_full[right] >= config.threshold_abs) right--;
    free(col_averages_full);
    if (left > right) { free(img->data); free(img); exit_error("Image seems blank horizontally.", NULL); }
    printf("  Horizontal content range: %d - %d\n", left, right);

    printf("Performing global vertical crop...\n");
    int global_top = 0; while (global_top < img->height && is_blank_row(img->data, img->width, global_top, img->height, left, right, config.threshold_abs)) global_top++;
    int global_bottom = img->height - 1; while (global_bottom >= global_top && is_blank_row(img->data, img->width, global_bottom, img->height, left, right, config.threshold_abs)) global_bottom--;
    if (global_top > global_bottom) { free(img->data); free(img); exit_error("Image seems blank vertically.", NULL); }
    int global_h = global_bottom - global_top + 1;
    printf("  Vertical content range: %d - %d (Height: %d)\n", global_top, global_bottom, global_h);

    printf("Finding track starts in ROI [%d-%d] (height %d) using initial smooth %d...\n", left, right, global_h, config.initial_smooth_window);
    int track_range_w = right - left + 1;
    if (track_range_w <= 0) { free(img->data); free(img); exit_error("ROI width is zero or negative.", NULL); }

    double *col_averages_roi = (double*)safe_malloc(track_range_w * sizeof(double));
    for (int c_rel = 0; c_rel < track_range_w; c_rel++) col_averages_roi[c_rel] = calculate_column_average_range(img->data, img->width, left + c_rel, global_top, global_bottom);
    double *smoothed_averages = (double*)safe_malloc(track_range_w * sizeof(double));
    smooth_1d_array(col_averages_roi, smoothed_averages, track_range_w, config.initial_smooth_window);
    free(col_averages_roi);

    int *track_starts = NULL; int starts_count = 0;
    int starts_capacity = 10;
    track_starts = (int*)safe_malloc(starts_capacity * sizeof(int));

    int first_col_is_dark = (track_range_w > 0 && smoothed_averages[0] < config.threshold_abs);
    if (first_col_is_dark) {
        track_starts[starts_count++] = 0;
        if (config.debug) printf("  Detected dark start at ROI beginning (offset 0 from col %d).\n", left);
    }
    int prev_is_dark = first_col_is_dark;
    for (int c_rel = 1; c_rel < track_range_w; c_rel++) {
        int is_dark = (smoothed_averages[c_rel] < config.threshold_abs);
        if (is_dark && !prev_is_dark) {
            if (starts_count >= starts_capacity) {
                starts_capacity *= 2;
                track_starts = (int*)safe_realloc(track_starts, starts_capacity * sizeof(int));
            }
            track_starts[starts_count++] = c_rel;
        }
        prev_is_dark = is_dark;
    }
    free(smoothed_averages);

    if (starts_count == 0) { free(track_starts); free(img->data); free(img); exit_error("No tracks found.", NULL); }
    printf("Found %d potential track start(s).\n", starts_count);
    if (config.debug) { printf("  Track start offsets relative to %d: ", left); for(int i=0; i<starts_count; ++i) printf("%d ", track_starts[i]); printf("\n"); }

    printf("Determining processing width and (V2) search delta...\n");
    config.calculated_avg_width = 0;
    if (starts_count > 1) {
        long long total_width_sum = 0; int segments_for_avg = 0;
        for (int i = 0; i < starts_count - 1; ++i) {
            int w = track_starts[i + 1] - track_starts[i];
            if (w > 0) { total_width_sum += w; segments_for_avg++; }
            else { fprintf(stderr,"Warning: Zero/negative distance between track starts %d and %d. Skipping for avg.\n", i+1, i+2); }
        }
        if (segments_for_avg > 0) config.calculated_avg_width = (int)round((double)total_width_sum / segments_for_avg);
    }
    printf("  Calculated average segment width: %d\n", config.calculated_avg_width);

    if (config.processing_width > 0) {
         config.final_processing_width = config.processing_width;
         printf("  Using fixed processing width (-w) as final: %d\n", config.final_processing_width);
    } else {
        int base_width = (config.calculated_avg_width > 0) ? config.calculated_avg_width : (track_range_w - track_starts[0]);
        base_width = (base_width <= 0) ? 1 : base_width;
        printf("  Using auto width base: %d\n", base_width);
        int reduced_width = (int)round((double)base_width * WIDTH_REDUCTION_FACTOR);
        config.final_processing_width = (reduced_width <= 0) ? 1 : reduced_width;
        printf("  Applying reduction -> final auto width: %d\n", config.final_processing_width);
    }
    if (config.final_processing_width <= 0) config.final_processing_width = 1;

    if (config.use_v2_logic) {
        if (config.search_delta <= 0) {
            config.search_delta = config.final_processing_width / DEFAULT_V2_SEARCH_DELTA_DIVIDER;
            if (config.search_delta <= 0) config.search_delta = 1;
            printf("  Using default V2 edge search delta (+/- %d)\n", config.search_delta);
        } else {
             printf("  Using specified V2 edge search delta (-d): +/- %d\n", config.search_delta);
        }
    }
    printf("  Final processing width to be used: %d\n", config.final_processing_width);
    printf("--------------------------\n");

    // ================================================================
    // Phase 2: Process Segments
    // ================================================================
    printf("\nPhase 2: Processing Segments (%s logic)...\n", config.use_v2_logic ? "V2 dynamic" : "V1 static");
    int segments_processed = 0;

    Image cropped_img_struct; // Buffer for V1 logic
    cropped_img_struct.data = NULL;
    if (!config.use_v2_logic) {
        cropped_img_struct.width = track_range_w;
        cropped_img_struct.height = global_h;
        cropped_img_struct.maxval = img->maxval;
        size_t cropped_size = (size_t)track_range_w * global_h;
        cropped_img_struct.data = (unsigned char*)safe_malloc(cropped_size);
        if (config.debug) printf("  (V1 mode: creating temporary cropped image buffer %dx%d)\n", track_range_w, global_h);
        for (int r = 0; r < global_h; r++) {
            memcpy(cropped_img_struct.data + (size_t)r * track_range_w,
                   img->data + (size_t)(global_top + r) * img->width + left,
                   track_range_w);
        }
    }

    // Loop through detected track starts
    for (int i = 0; i < starts_count; ++i) {
         int track_num = i + 1;
         int segment_start_rel = track_starts[i];
         int segment_end_rel = (i == starts_count - 1) ? (track_range_w - 1) : (track_starts[i+1] - 1);

         if (segment_end_rel < segment_start_rel) {
              fprintf(stderr, "Warning: Invalid boundaries for track %d (start_rel=%d, end_rel=%d). Skipping.\n",
                      track_num, segment_start_rel, segment_end_rel);
              continue;
         }

         // Dispatch to appropriate processing function
         if (config.use_v2_logic) {
             process_segment_v2(
                 img, global_top, global_bottom,
                 segment_start_rel, left,
                 track_num, &config, input_basename
             );
         } else {
             process_segment_v1(
                 cropped_img_struct.data, cropped_img_struct.width, cropped_img_struct.height,
                 segment_start_rel, segment_end_rel,
                 cropped_img_struct.maxval, // Use maxval from cropped struct
                 &config, track_num,
                 config.calculated_avg_width,
                 input_basename
             );
         }
         segments_processed++;
    }

    // ================================================================
    // Phase 3: Finalization
    // ================================================================
    printf("\nPhase 3: Finalizing...\n");
    printf("Processed %d segments.\n", segments_processed);
    if (segments_processed != starts_count) {
         printf("Warning: %d potential start(s) might have been skipped.\n", starts_count - segments_processed);
    }
    printf("Output segments saved in: %s\n", config.output_m_directory);

    // --- Cleanup ---
    printf("Cleaning up resources...\n");
    free(track_starts);
    if (cropped_img_struct.data != NULL) {
        free(cropped_img_struct.data);
    }
    free(img->data);
    free(img);

    printf("Done.\n");
    return EXIT_SUCCESS;
}

