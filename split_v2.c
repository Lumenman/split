#define _POSIX_C_SOURCE 200809L // Для getopt, strdup, basename
#include <stdio.h>
#include <stdlib.h> // Добавлено для qsort, exit
#include <string.h>
#include <math.h>
#include <unistd.h> // для getopt
#include <errno.h>
#include <ctype.h>   // Для isprint, tolower
#include <sys/stat.h> // Для mkdir, stat
#include <libgen.h>   // Для basename

// Для mkdir в Windows (если не используется POSIX-совместимая среда типа MinGW)
#ifdef _WIN32
#include <direct.h>
#define MKDIR(path, mode) _mkdir(path)
#else
#define MKDIR(path, mode) mkdir(path, mode)
#endif

#define MAX_PATH_LEN 512
#define MAX_LINE_LEN 256
#define PGM_MAXVAL_SUPPORTED 255

// Значения по умолчанию
#define DEFAULT_INITIAL_SMOOTH_WINDOW 1
#define DEFAULT_ROW_SMOOTH_WINDOW 1
#define DEFAULT_PROCESSING_WIDTH 0
#define DEFAULT_SEARCH_DELTA_DIVIDER 2
#define DEFAULT_MIN_RUN_LENGTH 1
#define DEFAULT_MEDIAN_WINDOW 3 // Размер окна медианного фильтра по умолчанию (-M)
#define WIDTH_REDUCTION_FACTOR 0.95

// --- Структуры ---

/**
 * @brief Структура для хранения параметров конфигурации.
 */
typedef struct {
    char input_filename[MAX_PATH_LEN];
    char output_m_directory[MAX_PATH_LEN];
    int threshold_abs;
    int initial_smooth_window; // -S
    int row_smooth_window;     // -s
    int median_window;         // -M (Размер окна медианного фильтра)
    int processing_width;      // -w
    int search_delta;          // -d
    int min_run_length;        // -k
    int debug;                 // -b
    // Расчетные значения
    int calculated_avg_width;
    int final_processing_width;
} Config;

/**
 * @brief Структура для хранения данных изображения PGM.
 */
typedef struct {
    unsigned char *data;
    int width;
    int height;
    int maxval;
} Image;

// --- Вспомогательные функции ---

void exit_error(const char *message, const char *details) {
    fprintf(stderr, "ERROR: %s%s%s\n", message, (details ? ": " : ""), (details ? details : ""));
    exit(EXIT_FAILURE);
}

void* safe_malloc(size_t size) {
    if (size == 0) size = 1;
    void *ptr = malloc(size);
    if (!ptr) exit_error("Memory allocation failed", strerror(errno));
    return ptr;
}

void* safe_realloc(void* ptr, size_t size) {
    if (size == 0) size = 1;
    void *new_ptr = realloc(ptr, size);
    if (!new_ptr && size > 0) exit_error("Memory reallocation failed", strerror(errno));
    return new_ptr;
}

// Функция сравнения для qsort
int int_compare(const void *a, const void *b) {
    int int_a = *((int*)a);
    int int_b = *((int*)b);
    if (int_a == int_b) return 0;
    else if (int_a < int_b) return -1;
    else return 1;
}


// --- Функции работы с PGM --- (Без изменений)

Image* read_pgm_p5(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) exit_error("Cannot open input PGM file", filename);

    Image *img = (Image*)safe_malloc(sizeof(Image));
    img->data = NULL;
    char line[MAX_LINE_LEN];
    int fields_read;
    int read_maxval = -1;

    if (!fgets(line, sizeof(line), fp) || strncmp(line, "P5", 2) != 0) {
        fclose(fp); free(img);
        exit_error("Input file is not a P5 PGM", filename);
    }

    fields_read = 0;
    do {
        if (!fgets(line, sizeof(line), fp)) { fclose(fp); free(img);
        exit_error("Could not read dimensions/maxval", filename); }
        if (line[0] != '#') fields_read = sscanf(line, "%d %d", &img->width, &img->height);
    } while (line[0] == '#' || fields_read < 2);

    if (fields_read != 2 || img->width <= 0 || img->height <= 0) {
        fclose(fp); free(img);
        exit_error("Invalid dimensions", filename);
    }

    fields_read = 0;
    do {
        if (!fgets(line, sizeof(line), fp)) { fclose(fp); free(img);
        exit_error("Could not read maxval", filename); }
        if (line[0] != '#') fields_read = sscanf(line, "%d", &read_maxval);
    } while (line[0] == '#' || fields_read < 1);

    if (fields_read != 1 || read_maxval <= 0 || read_maxval > PGM_MAXVAL_SUPPORTED) {
        fclose(fp);
        free(img);
        char err_details[100];
        snprintf(err_details, sizeof(err_details), "Value %d. Must be > 0 and <= %d.", read_maxval, PGM_MAXVAL_SUPPORTED);
        exit_error("Invalid or unsupported maxval", err_details);
    }
    img->maxval = read_maxval;


    size_t size = (size_t)img->width * img->height;
    if (size == 0) { fclose(fp); free(img); exit_error("Image size is zero", filename); }
    img->data = (unsigned char*)safe_malloc(size);

    int c = fgetc(fp);
    if (c == EOF) { fclose(fp); free(img->data); free(img); exit_error("EOF before pixel data", filename); }
    if (!isspace(c)) { ungetc(c, fp); }

    if (fread(img->data, 1, size, fp) != size) {
        if(ferror(fp)) { fclose(fp); free(img->data); free(img); exit_error("Failed to read pixel data", strerror(errno)); }
        else { fclose(fp); free(img->data); free(img); exit_error("Insufficient pixel data", filename); }
    }

    fclose(fp);
    return img;
}

int write_pgm_p5(const char *filename, const unsigned char *data, int width, int height, int maxval) {
    if (width <= 0 || height <= 0 || !data) {
        fprintf(stderr, "Warning: Invalid dimensions or data for PGM write (%dx%d) to %s\n", width, height, filename);
        return -1;
    }
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Warning: Cannot create output PGM file %s (%s)\n", filename, strerror(errno));
        return -1;
    }
    fprintf(fp, "P5\n%d %d\n%d\n", width, height, maxval > PGM_MAXVAL_SUPPORTED ? PGM_MAXVAL_SUPPORTED : maxval);
    size_t size = (size_t)width * height;
    if (fwrite(data, 1, size, fp) != size) {
        fprintf(stderr, "Warning: Failed to write all pixel data to %s\n", filename);
        fclose(fp);
        return -1;
    }
    fclose(fp);
    return 0;
}


// --- Функции анализа изображений --- (Без изменений)

int is_blank_row(const Image *img, int row, int col_start, int col_end, int threshold) {
    if (row < 0 || row >= img->height) return 1;
    col_start = (col_start < 0) ? 0 : col_start;
    col_end = (col_end >= img->width) ? img->width - 1 : col_end;
    if (col_start > col_end) return 1;
    size_t row_offset = (size_t)row * img->width;
    for (int c = col_start; c <= col_end; c++) {
        if (img->data[row_offset + c] < threshold) {
            return 0;
        }
    }
    return 1;
}

double calculate_column_average_range(const Image *img, int col, int row_start, int row_end) {
    if (col < 0 || col >= img->width || row_start < 0 || row_end >= img->height || row_start > row_end) {
        return (double)img->maxval;
    }
    long long sum = 0;
    int count = row_end - row_start + 1;
    for (int r = row_start; r <= row_end; r++) {
        sum += img->data[(size_t)r * img->width + col];
    }
    return (count > 0) ? (double)sum / count : (double)img->maxval;
}

void smooth_1d_array(const double *data_in, double *data_out, int n, int window_size) {
    if (window_size <= 1 || n == 0) {
        if (n > 0) memcpy(data_out, data_in, n * sizeof(double));
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

void smooth_horizontal_row(const Image *img, int row, int col_start, int col_end, int window_size, double *out_smoothed_row) {
     if (row < 0 || row >= img->height) return;
     col_start = (col_start < 0) ? 0 : col_start;
     col_end = (col_end >= img->width) ? img->width - 1 : col_end;
     if (col_start > col_end) return;
     int out_len = col_end - col_start + 1;
     if (window_size <= 1) {
         size_t row_offset = (size_t)row * img->width;
         for(int i=0; i<out_len; ++i) {
             out_smoothed_row[i] = (double)img->data[row_offset + col_start + i];
         }
         return;
     }
     int half_window_floor = window_size / 2;
     int half_window_ceil = (window_size - 1) / 2;
     size_t row_offset = (size_t)row * img->width;
     for (int c_idx = 0; c_idx < out_len; ++c_idx) {
         int c_orig = col_start + c_idx;
         double sum = 0;
         int count = 0;
         int k_start_orig = c_orig - half_window_floor;
         int k_end_orig = c_orig + half_window_ceil;
         for (int k_orig = k_start_orig; k_orig <= k_end_orig; ++k_orig) {
             if (k_orig >= 0 && k_orig < img->width) {
                 sum += img->data[row_offset + k_orig];
                 count++;
             }
         }
         out_smoothed_row[c_idx] = (count > 0) ? sum / count : (double)img->data[row_offset + c_orig];
     }
}

int find_left_edge_in_row(const Image *img, int r_img, int predicted_col_orig, const Config *config, double *smoothed_row_buffer) {
    if (r_img < 0 || r_img >= img->height) return -1;

    int search_start_orig = predicted_col_orig - config->search_delta;
    int search_end_orig = predicted_col_orig + config->search_delta;
    search_start_orig = (search_start_orig < 0) ? 0 : search_start_orig;
    search_end_orig = (search_end_orig >= img->width) ? img->width - 1 : search_end_orig;
    if (search_start_orig > search_end_orig) return -1;
    int search_len = search_end_orig - search_start_orig + 1;

    smooth_horizontal_row(img, r_img, search_start_orig, search_end_orig, config->row_smooth_window, smoothed_row_buffer);

    int found_edge_orig = -1;
    double prev_smoothed_val;

    if (search_start_orig > 0) {
        double temp_buf[1];
        smooth_horizontal_row(img, r_img, search_start_orig - 1, search_start_orig - 1, config->row_smooth_window, temp_buf);
        prev_smoothed_val = temp_buf[0];
    } else {
        prev_smoothed_val = (double)img->maxval;
    }

    for (int i = 0; i < search_len; ++i) {
        double current_smoothed_val = smoothed_row_buffer[i];
        int current_col_orig = search_start_orig + i;

        if (current_smoothed_val < config->threshold_abs && prev_smoothed_val >= config->threshold_abs) {
            if (config->min_run_length > 1) {
                int dark_run_count = 0;
                size_t row_offset = (size_t)r_img * img->width;
                for (int k = 0; k < config->min_run_length; ++k) {
                    int check_col_orig = current_col_orig + k;
                    if (check_col_orig < img->width) {
                         if (img->data[row_offset + check_col_orig] < config->threshold_abs) {
                             dark_run_count++;
                         } else { break; }
                    } else { break; }
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


// --- Вспомогательные функции для main ---

void get_pgm_basename(const char* filepath, char* out_basename, size_t max_len) {
    if (!filepath || !out_basename || max_len == 0) return;
    char filepath_copy[MAX_PATH_LEN];
    strncpy(filepath_copy, filepath, sizeof(filepath_copy) - 1);
    filepath_copy[sizeof(filepath_copy) - 1] = '\0';
    char *base = basename(filepath_copy);
    strncpy(out_basename, base, max_len - 1);
    out_basename[max_len - 1] = '\0';
    char *dot = strrchr(out_basename, '.');
    if (dot != NULL && dot != out_basename) {
        if (strlen(dot) == 4 &&
            (tolower(dot[1]) == 'p') && (tolower(dot[2]) == 'g') && (tolower(dot[3]) == 'm'))
        {
            *dot = '\0';
        }
    }
}

// --- main ---
int main(int argc, char *argv[]) {
    Config config;
    // --- Установка значений по умолчанию ---
    config.input_filename[0] = '\0';
    config.output_m_directory[0] = '\0';
    config.threshold_abs = -1;
    config.initial_smooth_window = DEFAULT_INITIAL_SMOOTH_WINDOW;
    config.row_smooth_window = DEFAULT_ROW_SMOOTH_WINDOW;
    config.median_window = DEFAULT_MEDIAN_WINDOW; // Устанавливаем значение по умолчанию для -M
    config.processing_width = DEFAULT_PROCESSING_WIDTH;
    config.search_delta = -1;
    config.min_run_length = DEFAULT_MIN_RUN_LENGTH;
    config.debug = 0;
    config.calculated_avg_width = 0;
    config.final_processing_width = 0;

    // --- Парсинг аргументов командной строки ---
    int opt;
    // Добавлена опция 'M'
    while ((opt = getopt(argc, argv, "i:m:t:S:s:M:w:d:k:b")) != -1) {
        switch (opt) {
            case 'i': strncpy(config.input_filename, optarg, MAX_PATH_LEN - 1); config.input_filename[MAX_PATH_LEN - 1] = '\0'; break;
            case 'm': strncpy(config.output_m_directory, optarg, MAX_PATH_LEN - 1); config.output_m_directory[MAX_PATH_LEN - 1] = '\0'; break;
            case 't': config.threshold_abs = atoi(optarg);
                      if (config.threshold_abs < 0 || config.threshold_abs > PGM_MAXVAL_SUPPORTED) exit_error("Threshold (-t) must be between 0 and 255", optarg);
                      break;
            case 'S': config.initial_smooth_window = atoi(optarg); if (config.initial_smooth_window < 1) exit_error("Initial smooth window (-S) must be >= 1", optarg); break;
            case 's': config.row_smooth_window = atoi(optarg); if (config.row_smooth_window < 1) exit_error("Row smooth window (-s) must be >= 1", optarg); break;
            case 'M': config.median_window = atoi(optarg);
                      if (config.median_window < 3 || config.median_window % 2 == 0) exit_error("Median filter window (-M) must be an odd integer >= 3", optarg);
                      break;
            case 'w': config.processing_width = atoi(optarg); if (config.processing_width < 0) exit_error("Processing width (-w) cannot be negative", optarg); break;
            case 'd': config.search_delta = atoi(optarg); if (config.search_delta <= 0) exit_error("Search delta (-d) must be positive", optarg); break;
            case 'k': config.min_run_length = atoi(optarg); if (config.min_run_length < 1) exit_error("Min run length (-k) must be >= 1", optarg); break;
            case 'b': config.debug = 1; break;
            case '?':
                 // Добавлена 'M' в строку проверки опций с аргументами
                 if (strchr("imtsSswdkM", optopt)) { fprintf(stderr, "Option -%c requires an argument.\n", optopt); }
                 else if (isprint(optopt)) { fprintf(stderr, "Unknown option `-%c'.\n", optopt); }
                 else { fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt); }
                 // Fall through
            default:
                 fprintf(stderr, "\nUsage: %s -i <input.pgm> -m <pgm_dir> -t <abs_thresh> [options]\n", argv[0]);
                 fprintf(stderr, "Required:\n");
                 fprintf(stderr, "  -i <input.pgm>    : Input PGM file (P5 format, maxval <= %d)\n", PGM_MAXVAL_SUPPORTED);
                 fprintf(stderr, "  -m <pgm_dir>    : Output directory to save aligned PGM segments\n");
                 fprintf(stderr, "  -t <threshold>  : Absolute brightness threshold [0-%d]\n", PGM_MAXVAL_SUPPORTED);
                 fprintf(stderr, "Optional Processing:\n");
                 fprintf(stderr, "  -S <smooth_I>   : Initial column average smoothing window (>=1, default: %d)\n", DEFAULT_INITIAL_SMOOTH_WINDOW);
                 fprintf(stderr, "  -s <smooth_R>   : Row horizontal smoothing window for edge tracking (>=1, default: %d)\n", DEFAULT_ROW_SMOOTH_WINDOW);
                 fprintf(stderr, "  -M <median_W>   : Median filter window size for edge tracking (odd >=3, default: %d)\n", DEFAULT_MEDIAN_WINDOW); // Опция -M
                 fprintf(stderr, "  -w <width>      : Segment processing width (>=0, 0=auto ~avg*%.2f, default: %d)\n", WIDTH_REDUCTION_FACTOR, DEFAULT_PROCESSING_WIDTH);
                 fprintf(stderr, "  -d <delta>      : Edge search window delta (+/-, optional, default: final_width/%d)\n", DEFAULT_SEARCH_DELTA_DIVIDER);
                 fprintf(stderr, "  -k <runlen>     : Min dark pixels run length to confirm edge (>=1, default: %d)\n", DEFAULT_MIN_RUN_LENGTH);
                 fprintf(stderr, "  -b              : Enable debug messages\n");
                 return EXIT_FAILURE;
        }
    }

    // --- Проверка обязательных аргументов ---
    if (strlen(config.input_filename) == 0) exit_error("Input PGM filename (-i) is required", NULL);
    if (strlen(config.output_m_directory) == 0) exit_error("Output PGM directory (-m) is required", NULL);
    if (config.threshold_abs == -1) exit_error("Parameter -t <threshold> (absolute value 0-255) is required", NULL);
    if (optind < argc) { fprintf(stderr, "Error: Unexpected non-option arguments found: %s ...\n", argv[optind]); return EXIT_FAILURE; }

    // --- Извлечение базового имени файла для PGM ---
    char input_basename[MAX_PATH_LEN];
    get_pgm_basename(config.input_filename, input_basename, sizeof(input_basename));

    // --- Вывод эффективных параметров ---
    printf("--- Effective Parameters ---\n");
    printf("Input PGM (-i): %s (basename: %s)\n", config.input_filename, input_basename);
    printf("Output PGM Dir (-m): %s\n", config.output_m_directory);
    printf("Threshold (-t): %d\n", config.threshold_abs);
    printf("Initial Smooth (-S): %d\n", config.initial_smooth_window);
    printf("Row Smooth (-s): %d\n", config.row_smooth_window);
    printf("Median Filter (-M): %d\n", config.median_window); // Вывод -M
    printf("Processing Width (-w): %d (0=auto with ~%.0f%% reduction)\n", config.processing_width, (1.0-WIDTH_REDUCTION_FACTOR)*100.0);
    printf("Min Run Length (-k): %d\n", config.min_run_length);
    printf("Debug (-b): %s\n", config.debug ? "Enabled" : "Disabled");

    // --- Создание папки для PGM сегментов (-m) ---
    struct stat st = {0};
    if (stat(config.output_m_directory, &st) == -1) {
        printf("Creating directory for PGM segments: %s\n", config.output_m_directory);
        int mkdir_ret;
        #ifdef _WIN32
            mkdir_ret = MKDIR(config.output_m_directory, 0);
        #else
            mkdir_ret = MKDIR(config.output_m_directory, 0755);
        #endif
        if (mkdir_ret != 0) exit_error("Could not create directory specified by -m", strerror(errno));
    } else if (!S_ISDIR(st.st_mode)) exit_error("Path specified by -m exists but is not a directory", config.output_m_directory);

    // ================================================================
    // Этап 1: Начальный анализ (без изменений)
    // ================================================================
    printf("\nPhase 1: Initial Analysis...\n");
    printf("Reading input PGM: %s...\n", config.input_filename);
    Image *img = read_pgm_p5(config.input_filename);
    printf("Read image: %d x %d, maxval=%d\n", img->width, img->height, img->maxval);

    printf("Finding overall content boundaries (reference only)...\n");
    double *col_averages_full = (double*)safe_malloc(img->width * sizeof(double));
    for (int c = 0; c < img->width; c++) col_averages_full[c] = calculate_column_average_range(img, c, 0, img->height - 1);
    int left = 0; while (left < img->width && col_averages_full[left] >= config.threshold_abs) left++;
    int right = img->width - 1; while (right > left && col_averages_full[right] >= config.threshold_abs) right--;
    free(col_averages_full);
    if (left >= right) { free(img->data); free(img); exit_error("Image seems blank horizontally based on column averages.", NULL); }
    printf("  Reference horizontal content range: %d - %d\n", left, right);

    printf("Performing global vertical crop...\n");
    int global_top = 0; while (global_top < img->height && is_blank_row(img, global_top, left, right, config.threshold_abs)) global_top++;
    int global_bottom = img->height - 1; while (global_bottom > global_top && is_blank_row(img, global_bottom, left, right, config.threshold_abs)) global_bottom--;
    if (global_top > global_bottom) { free(img->data); free(img); exit_error("Image seems blank vertically after checking content range.", NULL); }
    int global_h = global_bottom - global_top + 1;
    printf("  Global vertical content range: %d - %d (Height: %d)\n", global_top, global_bottom, global_h);

    printf("Finding track starts within reference range [%d-%d]...\n", left, right);
    int track_range_w = right - left + 1;
    if (track_range_w <= 0) { free(img->data); free(img); exit_error("Reference horizontal range width is zero or negative.", NULL); }
    double *col_averages_roi = (double*)safe_malloc(track_range_w * sizeof(double));
    for (int c_rel = 0; c_rel < track_range_w; c_rel++) col_averages_roi[c_rel] = calculate_column_average_range(img, left + c_rel, global_top, global_bottom);
    double *smoothed_averages = (double*)safe_malloc(track_range_w * sizeof(double));
    smooth_1d_array(col_averages_roi, smoothed_averages, track_range_w, config.initial_smooth_window);
    free(col_averages_roi);

    int *track_starts = NULL; int starts_count = 0; int starts_capacity = 10;
    track_starts = (int*)safe_malloc(starts_capacity * sizeof(int));

    int first_col_is_dark = 0;
    if (track_range_w > 0 && smoothed_averages[0] < config.threshold_abs) {
        track_starts[starts_count++] = 0;
        first_col_is_dark = 1;
        if (config.debug) printf("  Detected dark start at the very beginning (offset 0).\n");
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

    if (starts_count == 0) { free(track_starts); free(img->data); free(img); exit_error("No tracks found within the reference range.", NULL); }
    printf("Found %d potential track start(s) using initial smooth window %d.\n", starts_count, config.initial_smooth_window);
    if (config.debug) { printf("  Track start offsets relative to %d: ", left); for(int i=0; i<starts_count; ++i) printf("%d ", track_starts[i]); printf("\n"); }

    printf("Determining processing width and search delta...\n");
    config.calculated_avg_width = 0;
    if (starts_count > 1) {
        long long total_width_sum = 0; int segments_for_avg = 0;
        for (int i = 0; i < starts_count - 1; ++i) {
            int w = track_starts[i + 1] - track_starts[i];
            if (w > 0) { total_width_sum += w; segments_for_avg++; }
            else { fprintf(stderr,"Warning: Zero or negative distance between track starts %d and %d. Skipping for average width calculation.\n", i+1, i+2); }
        }
        if (segments_for_avg > 0) config.calculated_avg_width = (int)round((double)total_width_sum / segments_for_avg);
    }
    printf("  Calculated average segment width (Arithmetic Mean): %d\n", config.calculated_avg_width);

    int base_width;
    if (config.processing_width > 0) {
        base_width = config.processing_width;
        printf("  Using fixed processing width (-w) as base: %d\n", base_width);
    } else {
        if (config.calculated_avg_width > 0) {
            base_width = config.calculated_avg_width;
            printf("  Using calculated average processing width as base: %d\n", base_width);
        } else {
            base_width = track_range_w - track_starts[0]; // Fallback
            base_width = (base_width <= 0) ? 1 : base_width;
            printf("  Using fallback processing width as base: %d\n", base_width);
        }
        if (config.processing_width == 0) {
             int reduced_width = (int)round((double)base_width * WIDTH_REDUCTION_FACTOR);
             if (reduced_width <= 0) reduced_width = 1;
             printf("  Applying %.0f%% reduction to automatically calculated base width: %d -> %d\n", (1.0-WIDTH_REDUCTION_FACTOR)*100.0, base_width, reduced_width);
             base_width = reduced_width;
        }
    }

    config.final_processing_width = (base_width <= 0) ? 1 : base_width;

    if (config.search_delta <= 0) {
        config.search_delta = config.final_processing_width / DEFAULT_SEARCH_DELTA_DIVIDER;
        if (config.search_delta <= 0) config.search_delta = 1;
        printf("  Using default edge search delta (+/- %d) based on final width %d\n", config.search_delta, config.final_processing_width);
    } else {
        printf("  Using specified edge search delta (-d): +/- %d\n", config.search_delta);
    }
    printf("  Final processing width to be used: %d\n", config.final_processing_width);
    printf("--------------------------\n");

    // --- 1.6 Подготовка буфера для Этапа 2 ---
    int max_search_len = (config.search_delta * 2) + 1;
    double *smoothed_row_buffer = (double*)safe_malloc(max_search_len * sizeof(double));
    int pgm_segments_processed = 0;

    // ================================================================
    // Этап 2: Построчная обработка и сохранение сегментов
    // ================================================================
    printf("\nPhase 2: Processing Segments and saving PGM files...\n");
    for (int i = 0; i < starts_count; ++i) {
        int track_index = i + 1;
        int predicted_start_col_orig = left + track_starts[i]; // Начальное предсказание для первого ряда сегмента
        if (config.debug) printf(" Processing Track %d (predicted start col: %d)...\n", track_index, predicted_start_col_orig);

        // --- 2.1 Вертикальная обрезка ---
        int segment_top_local = global_top;
        int segment_bottom_local = global_bottom;
        int vc_col_start = predicted_start_col_orig - config.search_delta;
        int vc_col_end = (i == starts_count - 1) ?
                         (right + config.search_delta) :
                         (left + track_starts[i+1] - 1 + config.search_delta);
        vc_col_start = (vc_col_start < 0) ? 0 : vc_col_start;
        vc_col_end = (vc_col_end >= img->width) ? img->width - 1 : vc_col_end;
        while (segment_top_local <= segment_bottom_local && is_blank_row(img, segment_top_local, vc_col_start, vc_col_end, config.threshold_abs)) segment_top_local++;
        while (segment_bottom_local >= segment_top_local && is_blank_row(img, segment_bottom_local, vc_col_start, vc_col_end, config.threshold_abs)) segment_bottom_local--;
        int segment_h_local = segment_bottom_local - segment_top_local + 1;

        if (segment_h_local <= 0) {
            fprintf(stderr, "Warning: Track %d skipped. Calculated height is %d after vertical crop (%d-%d).\n", track_index, segment_h_local, segment_top_local, segment_bottom_local);
            continue;
        }
        if (config.debug) printf("  Track %d: Vertically cropped to rows %d-%d (height %d)\n", track_index, segment_top_local, segment_bottom_local, segment_h_local);

        // --- 2.3 Подготовка к сохранению PGM и фильтрации ---
        unsigned char *pgm_segment_data = NULL;
        char pgm_filename[MAX_PATH_LEN * 2];
        pgm_segment_data = (unsigned char*)safe_malloc((size_t)config.final_processing_width * segment_h_local);
        memset(pgm_segment_data, img->maxval > PGM_MAXVAL_SUPPORTED ? PGM_MAXVAL_SUPPORTED : img->maxval, (size_t)config.final_processing_width * segment_h_local);
        int snprintf_ret = snprintf(pgm_filename, sizeof(pgm_filename), "%s/%s_%04d.pgm", config.output_m_directory, input_basename, track_index);
        if (snprintf_ret < 0 || snprintf_ret >= sizeof(pgm_filename)) {
            fprintf(stderr, "Warning: Failed to create filename for PGM track %d. PGM saving disabled for this track.\n", track_index);
            free(pgm_segment_data);
            pgm_segment_data = NULL;
        }

        // --- Инициализация для медианного фильтра ---
        int *edge_history = (int*)safe_malloc(config.median_window * sizeof(int));
        int *sort_buffer = (int*)safe_malloc(config.median_window * sizeof(int)); // Буфер для сортировки
        int history_count = 0; // Количество валидных значений в истории
        int history_idx = 0;   // Индекс для добавления следующего значения (циклический)
        int current_median_edge = predicted_start_col_orig; // Используем начальное предсказание до первой медианы

        // --- 2.4 Цикл по строкам сегмента ---
        for (int r_img = segment_top_local; r_img <= segment_bottom_local; ++r_img) {
            int r_rel = r_img - segment_top_local;

            // --- 2.4.1 Поиск "сырого" края в текущей строке ---
            // Предсказанием служит медиана из ПРЕДЫДУЩЕЙ строки
            int raw_edge = find_left_edge_in_row(img, r_img, current_median_edge, &config, smoothed_row_buffer);

            // --- 2.4.2 Обновление истории краев ---
            int edge_to_store;
            if (raw_edge != -1) {
                edge_to_store = raw_edge; // Используем найденное значение
            } else {
                // Сбой поиска: используем результат медианы предыдущей строки
                // (или начальное предсказание, если это первый ряд)
                edge_to_store = current_median_edge;
                if (config.debug > 1) printf("  Track %d, Row %d: Edge search failed, using previous median %d\n", track_index, r_img, current_median_edge);
            }
            edge_history[history_idx] = edge_to_store;
            history_idx = (history_idx + 1) % config.median_window; // Сдвигаем индекс циклически
            if (history_count < config.median_window) {
                history_count++; // Увеличиваем счетчик, пока буфер не заполнится
            }

            // --- 2.4.3 Вычисление медианы ---
            if (history_count > 0) {
                // Копируем валидные значения из истории в буфер для сортировки
                // (history_count определяет, сколько значений копировать)
                // Логика копирования из циклического буфера:
                int sort_idx = 0;
                int start_copy_idx = (history_count == config.median_window) ? history_idx : 0; // Откуда начинать копирование в циклическом буфере
                for(int k=0; k<history_count; ++k) {
                    int current_copy_idx = (start_copy_idx + k) % config.median_window;
                    sort_buffer[sort_idx++] = edge_history[current_copy_idx];
                }

                // Сортируем буфер
                qsort(sort_buffer, history_count, sizeof(int), int_compare);

                // Медиана - средний элемент
                current_median_edge = sort_buffer[history_count / 2];
                 if (config.debug > 1) printf("  Track %d, Row %d: Raw=%d, Stored=%d, Median(%d points)=%d\n", track_index, r_img, raw_edge, edge_to_store, history_count, current_median_edge);
            }
             // Если history_count == 0 (самая первая строка, где поиск неудачен),
             // то current_median_edge останется predicted_start_col_orig

            // --- 2.4.4 Копирование данных для PGM, используя МЕДИАНУ ---
            if (pgm_segment_data) {
                int win_start_col_orig = current_median_edge; // <--- Используем медианный край!
                int win_end_col_orig = win_start_col_orig + config.final_processing_width - 1;
                int win_start_clamped = (win_start_col_orig < 0) ? 0 : win_start_col_orig;
                int win_end_clamped = (win_end_col_orig >= img->width) ? img->width - 1 : win_end_col_orig;
                int actual_width_to_read = win_end_clamped - win_start_clamped + 1;

                unsigned char *dest_row_ptr = pgm_segment_data + (size_t)r_rel * config.final_processing_width;

                if (actual_width_to_read > 0) {
                    size_t src_row_offset = (size_t)r_img * img->width;
                    unsigned char *src_ptr = img->data + src_row_offset + win_start_clamped;
                    int copy_len = (actual_width_to_read < config.final_processing_width) ?
                                   actual_width_to_read : config.final_processing_width;
                    memcpy(dest_row_ptr, src_ptr, copy_len);
                }
            }
        } // --- Конец цикла по строкам ---

        // --- 2.5 Сохранение PGM файла и очистка истории ---
        if (pgm_segment_data) {
            if (write_pgm_p5(pgm_filename, pgm_segment_data, config.final_processing_width, segment_h_local, img->maxval) == 0) {
                if(config.debug) printf("  Saved aligned PGM segment: %s (%dx%d)\n", pgm_filename, config.final_processing_width, segment_h_local);
            }
            free(pgm_segment_data);
        }
        free(edge_history); // Освобождаем буфер истории для этого трека
        free(sort_buffer);  // Освобождаем буфер сортировки
        pgm_segments_processed++;

    } // --- Конец цикла по сегментам ---

    free(smoothed_row_buffer);

    // ================================================================
    // Этап 3: Финализация (только очистка)
    // ================================================================
    printf("\nPhase 3: Finalizing...\n");
    printf("Processed %d valid segments.\n", pgm_segments_processed);
    if (pgm_segments_processed != starts_count) {
         printf("Warning: %d segment(s) were skipped due to zero height after vertical crop.\n", starts_count - pgm_segments_processed);
    }
    printf("Aligned PGM segments saved in: %s\n", config.output_m_directory);

    // --- Очистка ---
    printf("Cleaning up resources...\n");
    free(track_starts);
    free(img->data);
    free(img);

    printf("Done.\n");
    return EXIT_SUCCESS;
}
