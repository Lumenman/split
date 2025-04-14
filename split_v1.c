#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>
#include <dirent.h> // Для работы с директориями (POSIX)
#include <unistd.h> // Для getopt
#include <ctype.h>  // Для isprint
#include <math.h>   // Для round

#define MAX_PATH_LEN 512
#define MAX_LINE_LEN 256
#define PGM_MAXVAL_SUPPORTED 255 // Максимальное значение яркости, поддерживаемое этим кодом
#define WIDTH_REDUCTION_FACTOR 0.95 // Коэффициент уменьшения авто-ширины

// --- Структуры ---

/**
 * @brief Структура для хранения параметров конфигурации.
 */
typedef struct {
    char input_filename[MAX_PATH_LEN]; ///< Путь к входному PGM файлу (-i). Обязательный.
    char output_folder[MAX_PATH_LEN];  ///< Папка для сохранения результатов (-m). Обязательная.
    int threshold_abs;                 ///< Абсолютный порог яркости (-t, 0-255). Обязательный.
    int smooth_window;                 ///< Размер окна сглаживания (-s).
    int track_width;                   ///< Фиксированная ширина трека (-w). 0 = авто.
} Config;

/**
 * @brief Структура для хранения данных изображения PGM.
 */
typedef struct {
    unsigned char *data;  ///< Массив пикселей изображения (8-бит grayscale).
    int width;            ///< Ширина изображения.
    int height;           ///< Высота изображения.
    int maxval;           ///< Максимальное значение яркости (прочитанное из файла).
} Image;

// --- Вспомогательные функции ---

/**
 * @brief Выводит сообщение об ошибке в stderr и завершает программу.
 * @param message Основное сообщение об ошибке.
 * @param details Дополнительные детали (например, strerror(errno) или имя файла).
 */
void exit_error(const char *message, const char *details) {
    fprintf(stderr, "ERROR: %s%s%s\n",
            message ? message : "",
            details ? " (" : "",
            details ? details : ")");
    exit(EXIT_FAILURE);
}

/**
 * @brief Безопасное выделение памяти с проверкой результата.
 * Завершает программу в случае ошибки выделения.
 * @param size Размер выделяемой памяти в байтах.
 * @return Указатель на выделенную память.
 */
void* safe_malloc(size_t size) {
    if (size == 0) size = 1;
    void *ptr = malloc(size);
    if (!ptr) {
        exit_error("Memory allocation failed", strerror(errno));
    }
    return ptr;
}

/**
 * @brief Безопасное перераспределение памяти с проверкой результата.
 * Завершает программу в случае ошибки перераспределения.
 * @param ptr Указатель на ранее выделенную память (или NULL).
 * @param size Новый размер выделяемой памяти в байтах.
 * @return Указатель на перераспределенную память.
 */
void* safe_realloc(void* ptr, size_t size) {
    if (size == 0) size = 1;
    void *new_ptr = realloc(ptr, size);
    if (!new_ptr && size > 0) {
        exit_error("Memory reallocation failed", strerror(errno));
    }
    return new_ptr;
}


// --- Функции чтения/записи PGM ---

/**
 * @brief Читает PGM-файл (P5 - бинарный) с проверками безопасности.
 * Поддерживает только maxval <= PGM_MAXVAL_SUPPORTED.
 * @param filename Путь к входному PGM файлу.
 * @return Указатель на структуру Image с данными или вызывает exit_error при ошибке.
 */
Image* read_pgm_p5(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        exit_error("Cannot open input PGM file", filename);
    }

    Image *img = (Image*)safe_malloc(sizeof(Image));
    img->data = NULL;
    char line[MAX_LINE_LEN];
    int fields_read;

    // Проверка формата P5
    if (!fgets(line, sizeof(line), fp) || strncmp(line, "P5", 2) != 0) {
        fclose(fp);
        free(img);
        exit_error("Input file is not a P5 PGM", filename);
    }

    // Чтение размеров, пропуская комментарии
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
        fclose(fp);
        free(img);
        exit_error("Invalid dimensions in PGM file", filename);
    }

    // Чтение maxval, пропуская комментарии
    fields_read = 0;
    do {
        if (!fgets(line, sizeof(line), fp)) {
            fclose(fp); free(img);
            exit_error("Could not read maxval from PGM header", filename);
        }
        if (line[0] != '#') {
            fields_read = sscanf(line, "%d", &img->maxval);
        }
    } while (line[0] == '#' || fields_read < 1);

    if (fields_read != 1 || img->maxval <= 0 ) {
        fclose(fp);
        free(img);
        exit_error("Invalid maxval in PGM file", filename);
    }
    if (img->maxval > PGM_MAXVAL_SUPPORTED) {
        fprintf(stderr, "Warning: PGM maxval (%d) in %s exceeds supported value (%d).\n",
                img->maxval, filename, PGM_MAXVAL_SUPPORTED);
        fprintf(stderr, "         Pixel data will be read, but results might be unexpected.\n");
    }

    // Чтение данных изображения
    size_t size = (size_t)img->width * img->height;
    if (size == 0) {
        fclose(fp); free(img);
        exit_error("Image size is zero based on dimensions", filename);
    }
    img->data = (unsigned char*)safe_malloc(size);

    int c = fgetc(fp);
    if (c == EOF) {
        fclose(fp); free(img->data); free(img);
        exit_error("Unexpected EOF before pixel data", filename);
    }
    if (!isspace(c)) {
         ungetc(c, fp);
    }

    if (fread(img->data, 1, size, fp) != size) {
        if (ferror(fp)) {
            fclose(fp); free(img->data); free(img);
            exit_error("Failed to read pixel data", strerror(errno));
        } else {
            fclose(fp); free(img->data); free(img);
            exit_error("Insufficient pixel data in file", filename);
        }
    }

    fclose(fp);
    return img;
}

/**
 * @brief Записывает PGM-файл (P5 - бинарный) с проверками.
 * @param filename Путь к выходному файлу.
 * @param img Указатель на структуру Image с данными для записи.
 * @return 0 при успехе, -1 при ошибке (выводится предупреждение).
 */
int write_pgm_p5(const char *filename, const Image *img) {
    if (!img || !img->data || img->width <= 0 || img->height <= 0) {
         fprintf(stderr, "Warning: Invalid image data provided for writing to %s.\n", filename);
         return -1;
    }

    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Warning: Cannot create output PGM file %s (%s).\n", filename, strerror(errno));
        return -1;
    }

    fprintf(fp, "P5\n%d %d\n%d\n", img->width, img->height, img->maxval > PGM_MAXVAL_SUPPORTED ? PGM_MAXVAL_SUPPORTED : img->maxval);

    size_t size = (size_t)img->width * img->height;
    if (fwrite(img->data, 1, size, fp) != size) {
        fprintf(stderr, "Warning: Failed to write all pixel data to %s.\n", filename);
        fclose(fp);
        return -1;
    }

    fclose(fp);
    return 0;
}


// --- Функции обработки изображений ---

/**
 * @brief Вычисляет среднее значение яркости в столбце исходного изображения.
 * @param img Указатель на структуру Image.
 * @param col Индекс столбца.
 * @return Среднее значение яркости или значение img->maxval при некорректном столбце.
 */
double calculate_column_average(const Image *img, int col) {
    if (col < 0 || col >= img->width || img->height <= 0) {
        return (double)img->maxval;
    }
    long long sum = 0;
    for (int i = 0; i < img->height; i++) {
        sum += img->data[(size_t)i * img->width + col];
    }
    return (img->height > 0) ? (double)sum / img->height : (double)img->maxval;
}

/**
 * @brief Проверяет наличие темного пикселя (яркость < threshold) в столбце сегмента.
 * @param track_data Указатель на данные сегмента (часть основного изображения или копия).
 * @param track_w Ширина сегмента.
 * @param track_h Высота сегмента.
 * @param col Индекс столбца внутри сегмента.
 * @param threshold Абсолютный порог яркости.
 * @return 1 если темный пиксель найден, 0 в противном случае.
 */
int column_has_dark_pixel(const unsigned char* track_data, int track_w, int track_h, int col, int threshold) {
    if (col < 0 || col >= track_w || track_h <= 0) {
        return 0;
    }
    for (int i = 0; i < track_h; i++) {
        if (track_data[(size_t)i * track_w + col] < threshold) {
            return 1;
        }
    }
    return 0;
}

/**
 * @brief Проверяет, является ли строка изображения "пустой" (не содержит темных пикселей)
 * в заданном диапазоне столбцов.
 * @param img_data Указатель на данные всего изображения.
 * @param img_width Ширина всего изображения.
 * @param row Индекс проверяемой строки.
 * @param col_start Начальный столбец диапазона (включительно).
 * @param col_end Конечный столбец диапазона (включительно).
 * @param threshold Абсолютный порог яркости.
 * @return 1 если строка пустая в диапазоне, 0 если содержит темный пиксель.
 */
int is_blank_row(const unsigned char *img_data, int img_width, int row, int col_start, int col_end, int threshold) {
    if (col_start > col_end || col_start >= img_width || col_end < 0) return 1;
    col_start = (col_start < 0) ? 0 : col_start;
    col_end = (col_end >= img_width) ? img_width - 1 : col_end;

    size_t row_offset = (size_t)row * img_width;

    for (int i = col_start; i <= col_end; i++) {
        if (img_data[row_offset + i] < threshold) {
            return 0;
        }
    }
    return 1;
}

/**
 * @brief Сглаживает 1D массив данных скользящим средним.
 * @param data_in Входной массив данных.
 * @param data_out Выходной массив для сглаженных данных (должен быть того же размера).
 * @param n Размер массивов.
 * @param window_size Размер окна сглаживания (должен быть >= 1).
 */
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


// --- Функция обработки и сохранения СЕГМЕНТА ---

/**
 * @brief Выделяет сегмент изображения, обрезает пустые края, применяет фиксированную
 * или автоматическую ширину (с фактором 0.95) и сохраняет результат в PGM-файл.
 * @param cropped_data Указатель на пиксели предварительно обрезанного общего изображения.
 * @param cropped_w Ширина предварительно обрезанного изображения.
 * @param cropped_h Высота предварительно обрезанного изображения.
 * @param segment_start_col Начальный столбец сегмента (в координатах cropped_data).
 * @param segment_end_col Конечный столбец сегмента (в координатах cropped_data).
 * @param original_maxval Исходное значение maxval (для заголовка PGM).
 * @param config Указатель на структуру конфигурации (порог, ширина и т.д.).
 * @param track_num Порядковый номер трека (для имени файла).
 * @param average_width_hint Рассчитанная средняя ширина сегмента (используется, если -w=0).
 */
void extract_save_segment(
    const unsigned char* cropped_data, int cropped_w, int cropped_h,
    int segment_start_col, int segment_end_col,
    int original_maxval, const Config* config,
    int track_num,
    int average_width_hint
) {
    // 1. Проверка и выделение данных исходного сегмента (копирование)
    int initial_segment_w = segment_end_col - segment_start_col + 1;
    if (initial_segment_w <= 0 || cropped_h <= 0) {
        fprintf(stderr,"Warning: Invalid initial segment dimensions for track %d (%dx%d). Skipping.\n",
                track_num, initial_segment_w, cropped_h);
        return;
    }

    Image segment_img;
    segment_img.width = initial_segment_w;
    segment_img.height = cropped_h;
    segment_img.maxval = original_maxval;
    segment_img.data = (unsigned char*)safe_malloc((size_t)segment_img.width * segment_img.height);

    for (int r = 0; r < cropped_h; r++) {
        memcpy(segment_img.data + (size_t)r * segment_img.width,
               cropped_data + (size_t)r * cropped_w + segment_start_col,
               segment_img.width);
    }

    // 2. Уточнение ЛЕВОГО края
    int left_crop_col_idx = 0;
    while(left_crop_col_idx < segment_img.width) {
        if (column_has_dark_pixel(segment_img.data, segment_img.width, segment_img.height, left_crop_col_idx, config->threshold_abs)) {
            break;
        }
        left_crop_col_idx++;
    }

    int right_boundary_in_segment = segment_img.width - 1;
    if (left_crop_col_idx > right_boundary_in_segment) {
        if (config->smooth_window > 1) {
             printf("Skipped potentially empty track %d after left-crop check (col %d).\n", track_num, segment_start_col + left_crop_col_idx);
        }
        free(segment_img.data); return;
    }

    int content_w_after_left_crop = right_boundary_in_segment - left_crop_col_idx + 1;
    if (content_w_after_left_crop <= 0) {
        printf("Skipped track %d due to zero width after left refinement.\n", track_num);
        free(segment_img.data); return;
    }

    // 3. Вертикальная обрезка
    int t_top = 0;
    int t_bottom = segment_img.height - 1;
    while (t_top < segment_img.height && is_blank_row(segment_img.data, segment_img.width, t_top, left_crop_col_idx, right_boundary_in_segment, config->threshold_abs)) {
        t_top++;
    }
    while (t_bottom > t_top && is_blank_row(segment_img.data, segment_img.width, t_bottom, left_crop_col_idx, right_boundary_in_segment, config->threshold_abs)) {
        t_bottom--;
    }

    int final_h = t_bottom - t_top + 1;
    if (final_h <= 0) {
        printf("Skipped track %d after vertical crop (rows %d-%d).\n", track_num, t_top, t_bottom);
        free(segment_img.data);
        return;
    }

    // 4. Определение финальной ширины холста для сохранения
    int final_w;
    int actual_content_width = content_w_after_left_crop;
    char width_source[50]; // Для информационного сообщения

    if (config->track_width > 0) {
        // Приоритет у -w
        final_w = config->track_width;
        snprintf(width_source, sizeof(width_source), "forced by -w");
    } else if (average_width_hint > 0) {
        // Используем среднюю ширину, умноженную на фактор
        int reduced_width = (int)round((double)average_width_hint * WIDTH_REDUCTION_FACTOR);
        if (reduced_width <= 0) reduced_width = 1; // Ширина должна быть минимум 1
        final_w = reduced_width;
        snprintf(width_source, sizeof(width_source), "avg*%.2f", WIDTH_REDUCTION_FACTOR);
    } else {
        // Используем текущую ширину контента как fallback
        final_w = actual_content_width;
        snprintf(width_source, sizeof(width_source), "content width");
    }

    if (final_w <= 0) {
        printf("Warning: Final width is zero or negative for track %d. Skipping.\n", track_num);
        free(segment_img.data); return;
    }

    // 5. Создание финального изображения и копирование данных
    Image final_track_img;
    final_track_img.width = final_w;
    final_track_img.height = final_h;
    final_track_img.maxval = original_maxval;
    final_track_img.data = (unsigned char*)safe_malloc((size_t)final_track_img.width * final_track_img.height);

    memset(final_track_img.data,
           (original_maxval > PGM_MAXVAL_SUPPORTED ? PGM_MAXVAL_SUPPORTED : original_maxval),
           (size_t)final_track_img.width * final_track_img.height);

    int copy_w = (final_w < actual_content_width) ? final_w : actual_content_width;

    if (copy_w > 0) {
        for (int r = 0; r < final_h; r++) {
            unsigned char* source_ptr = segment_img.data + (size_t)(t_top + r) * segment_img.width + left_crop_col_idx;
            unsigned char* dest_ptr = final_track_img.data + (size_t)r * final_track_img.width;
            memcpy(dest_ptr, source_ptr, copy_w);
        }
    }

    // 6. Сохранение отдельного файла трека
    char fname[MAX_PATH_LEN * 2];
    int chars_written = snprintf(fname, sizeof(fname), "%s/track_%03d.pgm", config->output_folder, track_num);

    if (chars_written < 0 || chars_written >= sizeof(fname)) {
         fprintf(stderr, "Warning: Failed to format filename for track %d. Skipping save.\n", track_num);
    } else {
        if (write_pgm_p5(fname, &final_track_img) == 0) {
             // Используем сохраненный источник ширины для сообщения
             printf("Saved track: %s (width %d - %s)\n", fname, final_track_img.width, width_source);
        }
    }

    // 7. Освобождение памяти
    free(segment_img.data);
    free(final_track_img.data);
}


int main(int argc, char *argv[]) {
    Config config;

    // --- Установка значений по умолчанию ---
    config.input_filename[0] = '\0'; // -i обязателен
    config.output_folder[0] = '\0';  // -m обязателен
    config.threshold_abs = -1;       // -t обязателен
    config.smooth_window = 1;
    config.track_width = 0;

    // --- Парсинг опций командной строки ---
    int opt;
    // Добавлена опция 'i'
    while ((opt = getopt(argc, argv, "i:m:s:t:w:")) != -1) {
         switch (opt) {
             case 'i': // Входной файл
                 strncpy(config.input_filename, optarg, sizeof(config.input_filename) - 1);
                 config.input_filename[sizeof(config.input_filename) - 1] = '\0';
                 break;
             case 'm': // Папка вывода
                 strncpy(config.output_folder, optarg, sizeof(config.output_folder) - 1);
                 config.output_folder[sizeof(config.output_folder) - 1] = '\0';
                 break;
             case 's': // Окно сглаживания
                 config.smooth_window = atoi(optarg);
                 if (config.smooth_window <= 0) {
                     exit_error("Smooth window (-s) must be a positive integer", optarg);
                 }
                 break;
             case 't': // Абсолютный порог
                 config.threshold_abs = atoi(optarg);
                 if (config.threshold_abs < 0 || config.threshold_abs > PGM_MAXVAL_SUPPORTED) {
                     char err_details[100];
                     snprintf(err_details, sizeof(err_details), "Value '%s'. Must be between 0 and %d.", optarg, PGM_MAXVAL_SUPPORTED);
                     exit_error("Invalid threshold (-t)", err_details);
                 }
                 break;
             case 'w': // Фиксированная ширина
                 config.track_width = atoi(optarg);
                 if (config.track_width < 0) {
                     exit_error("Track width (-w) cannot be negative", optarg);
                 }
                 break;
             case '?': // Обработка ошибок getopt
                 // Обновляем сообщение об использовании
                 fprintf(stderr, "\nUsage: %s -i <input.pgm> -m <output_dir> -t <threshold> [-s smooth] [-w width]\n", argv[0]);
                 fprintf(stderr, "  -i <input.pgm>     : Input PGM file (P5 format, maxval <= %d, required)\n", PGM_MAXVAL_SUPPORTED);
                 fprintf(stderr, "  -m <output_dir>  : Directory to save output tracks (required)\n");
                 fprintf(stderr, "  -t <threshold>   : Absolute brightness threshold [0-%d] (required)\n", PGM_MAXVAL_SUPPORTED);
                 fprintf(stderr, "  -s <smooth>      : Smoothing window size (positive integer, default=1)\n");
                 fprintf(stderr, "  -w <width>       : Force output track width (>=0, 0=auto with ~%.0f%% reduction, default=0)\n", (1.0-WIDTH_REDUCTION_FACTOR)*100.0);
                 return EXIT_FAILURE;
             default:
                 abort();
         }
     }

    // --- Проверка обязательных аргументов ---
    if (strlen(config.input_filename) == 0) {
        exit_error("Input PGM filename (-i) is required.", NULL);
    }
    if (strlen(config.output_folder) == 0) {
        exit_error("Output folder (-m) is required.", NULL);
    }
    if (config.threshold_abs == -1) {
        exit_error("Absolute threshold (-t) is required.", NULL);
    }

    // Проверка наличия лишних аргументов не-опций
    if (optind < argc) {
        exit_error("Unexpected non-option arguments found.", argv[optind]);
    }


    // --- Вывод используемых параметров ---
    printf("--- Input Parameters ---\n");
    printf("Input scan file (-i): %s\n", config.input_filename); // Используем config.input_filename
    printf("Output folder (-m)  : %s\n", config.output_folder);
    printf("Smooth window (-s)  : %d\n", config.smooth_window);
    printf("Threshold (-t)      : %d\n", config.threshold_abs);
    printf("Forced width (-w)   : %d (0 = auto width ~avg*%.2f)\n", config.track_width, WIDTH_REDUCTION_FACTOR); // Обновлено описание -w
    printf("------------------------\n");

    // --- Создание выходной папки ---
    struct stat st = {0};
    if (stat(config.output_folder, &st) == -1) {
        if (mkdir(config.output_folder, 0755) != 0) {
            exit_error("Could not create output folder", config.output_folder);
        }
        printf("Created output directory: %s\n", config.output_folder);
    } else if (!S_ISDIR(st.st_mode)) {
        exit_error("Output path exists but is not a directory", config.output_folder);
    }

    // --- Чтение исходного изображения ---
    printf("Reading input image...\n");
    Image *img = read_pgm_p5(config.input_filename); // Используем config.input_filename
    printf("Read main image: %s (%d x %d, maxval=%d)\n",
           config.input_filename, img->width, img->height, img->maxval); // Используем config.input_filename

    // --- Вычисление средних яркостей столбцов ---
    if (img->width <= 0) exit_error("Image width is zero or negative after reading.", NULL);
    double *col_averages = (double*)safe_malloc(img->width * sizeof(double));
    for(int c = 0; c < img->width; c++) {
        col_averages[c] = calculate_column_average(img, c);
    }

    // --- Первичная обрезка пустых краев всего изображения ---
    printf("Cropping blank borders...\n");
    int left = 0;
    int right = img->width - 1;
    while (left < img->width && col_averages[left] >= config.threshold_abs) {
        left++;
    }
    while (right > left && col_averages[right] >= config.threshold_abs) {
        right--;
    }
    if (left > right) {
        printf("Image seems blank based on column averages and threshold %d.\n", config.threshold_abs);
        free(col_averages); free(img->data); free(img);
        return EXIT_SUCCESS;
    }

    int top = 0;
    int bottom = img->height - 1;
    while (top < img->height && is_blank_row(img->data, img->width, top, left, right, config.threshold_abs)) {
        top++;
    }
    while (bottom > top && is_blank_row(img->data, img->width, bottom, left, right, config.threshold_abs)) {
        bottom--;
    }
    if (top > bottom) {
        printf("Image seems blank based on rows after column crop (threshold %d).\n", config.threshold_abs);
        free(col_averages); free(img->data); free(img);
        return EXIT_SUCCESS;
    }

    int cropped_w = right - left + 1;
    int cropped_h = bottom - top + 1;
    if (cropped_w <= 0 || cropped_h <= 0) {
        exit_error("Cropped dimensions are zero or negative after border removal.", NULL);
    }
    printf("Cropped image dimensions: %d x %d (Cols %d-%d, Rows %d-%d)\n",
           cropped_w, cropped_h, left, right, top, bottom);

    // --- Создание буфера для обрезанной части изображения ---
    Image cropped_img_struct;
    cropped_img_struct.width = cropped_w;
    cropped_img_struct.height = cropped_h;
    cropped_img_struct.maxval = img->maxval;
    size_t cropped_size = (size_t)cropped_w * cropped_h;
    cropped_img_struct.data = (unsigned char*)safe_malloc(cropped_size);

    for (int r = 0; r < cropped_h; r++) {
        memcpy(cropped_img_struct.data + (size_t)r * cropped_w,
               img->data + (size_t)(top + r) * img->width + left,
               cropped_w);
    }
    int original_maxval = img->maxval;
    free(img->data); free(img);
    // col_averages все еще нужны!

    // --- Сглаживание средних яркостей (обрезанной области) ---
    if (cropped_w <= 0) exit_error("Cropped width is zero or negative before smoothing.", NULL);
    double *smooth_avg = (double*)safe_malloc(cropped_w * sizeof(double));
    printf("Smoothing column averages (window size %d)...\n", config.smooth_window);

    double *cropped_col_averages = (double*)safe_malloc(cropped_w * sizeof(double));
    for(int c = 0; c < cropped_w; ++c) {
        cropped_col_averages[c] = col_averages[left + c];
    }
    smooth_1d_array(cropped_col_averages, smooth_avg, cropped_w, config.smooth_window);

    free(cropped_col_averages);
    free(col_averages);

    // --- Поиск начал всех треков ---
    printf("Finding track starts...\n");
    int *track_starts = NULL;
    int starts_count = 0;
    int starts_capacity = 10;
    track_starts = (int*)safe_malloc(starts_capacity * sizeof(int));

    for (int col = 0; col < cropped_w; col++) {
        int is_dark = smooth_avg[col] < config.threshold_abs;
        int prev_is_dark = (col > 0) ? (smooth_avg[col - 1] < config.threshold_abs) : !is_dark;

        if (is_dark && !prev_is_dark) {
            if (starts_count >= starts_capacity) {
                starts_capacity *= 2;
                track_starts = (int*)safe_realloc(track_starts, starts_capacity * sizeof(int));
            }
            track_starts[starts_count++] = col;
        }
    }
    free(smooth_avg);

    if (starts_count == 0) {
        printf("No tracks found based on threshold %d and smoothing %d.\n",
               config.threshold_abs, config.smooth_window);
        free(cropped_img_struct.data); free(track_starts);
        return EXIT_SUCCESS;
    }
    printf("Found %d potential track start(s).\n", starts_count);

    // --- Вычисление средней ширины сегментов (для опции -w 0) ---
    int average_segment_width = 0;
    if (starts_count > 1) {
        long long total_segment_width = 0;
        int num_segments_for_avg = 0;
        for (int i = 0; i < starts_count - 1; ++i) {
             int segment_w = track_starts[i+1] - track_starts[i];
             if (segment_w > 0) {
                 total_segment_width += segment_w;
                 num_segments_for_avg++;
             } else {
                  fprintf(stderr, "Warning: Inconsistent track start order or zero width segment (%d -> %d at col %d). Skipping segment in average calculation.\n",
                          i+1, i+2, track_starts[i]);
             }
        }
        if (num_segments_for_avg > 0) {
             average_segment_width = (int)round((double)total_segment_width / num_segments_for_avg);
             if (average_segment_width <= 0) average_segment_width = 1;
             printf("Average segment width (basis for -w 0): %d\n", average_segment_width);
        } else {
             printf("Could not calculate average width (no valid segments between tracks).\n");
        }
    } else {
        printf("Only one track found, average width not applicable.\n");
    }

    // --- Обработка и сохранение каждого сегмента ---
    printf("Processing segments and saving individual tracks...\n");
    int num_tracks_saved = 0;
    for (int i = 0; i < starts_count; ++i) {
        int segment_start = track_starts[i];
        int segment_end = (i == starts_count - 1) ? (cropped_w - 1) : (track_starts[i+1] - 1);

        if (segment_end < segment_start) {
             fprintf(stderr, "Warning: Invalid segment boundaries for track %d (start=%d, end=%d). Skipping.\n",
                     i + 1, segment_start, segment_end);
             continue;
        }

        extract_save_segment(
            cropped_img_struct.data, cropped_img_struct.width, cropped_img_struct.height,
            segment_start, segment_end, original_maxval, &config, i + 1,
            average_segment_width // Передаем СРЕДНЮЮ ширину как подсказку
        );
        num_tracks_saved++;
    }

    // --- Финальная очистка ---
    free(cropped_img_struct.data);
    free(track_starts);

    printf("Finished. Processed %d segments.\n", num_tracks_saved);
    printf("Individual tracks saved in: %s\n", config.output_folder);
    return EXIT_SUCCESS;
}
