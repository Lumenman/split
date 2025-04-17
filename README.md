# Split PGM into segments

Этот проект предоставляет инструменты для разделения PGM-изображений на сегменты (дорожки). Он связан с другим [проектом](http://zenpho.co.uk/paper.shtml), который позволяет "распечатать" звук на бумагу в форме оптической звуковой дорожки, состоящей из ряда сегментов на листе бумаги А4.

## Функциональность

Для разделения исходного изображения на сегменты используется программа, скомпилированная из исходного кода `split.c` (или как будет назван исполняемый файл).

Эта программа объединяет **две основные логики обработки (V1 и V2)**, которые можно выбрать с помощью параметра командной строки `-e`. По умолчанию используется логика V1.

### Режим V1 (по умолчанию, без флага `-e`)

* **Начальное обнаружение:** Находит примерные границы дорожек по средним значениям яркости столбцов исходного изображения. Начальное сглаживание этих средних значений управляется опцией `-s`.
* **Определение статического края:** Для каждого найденного сегмента определяет его левую границу *один раз*, находя первый столбец, содержащий темный пиксель (яркость ниже порога `-t`).
* **Обрезка:** Обрезает общие пустые поля сверху и снизу изображения, а затем каждый сегмент обрезается по вертикали на основе найденного *статического* левого края.
* **Сохранение:** Копирует данные сегмента, начиная от найденного статического левого края, в выходной файл. Ширина выходного файла задается опцией `-w` или рассчитывается автоматически.
* **Ограничение:** Этот алгоритм не "следит" за левым краем вдоль всей дорожки. При наклоне или кривизне исходного скана это может привести к обрезке полезной части информации слева в некоторых частях дорожки.

### Режим V2 (активируется флагом `-e`)

* **Начальное обнаружение:** Аналогично V1, находит примерные границы дорожек. Начальное сглаживание управляется опцией `-S`.
* **Динамическое отслеживание края:** Этот режим *построчно* отслеживает левый край каждой дорожки:
    * Для каждой строки предсказывается положение края (на основе предыдущей строки).
    * В небольшом окне (`-d`) вокруг предсказания ищется фактический переход от светлого к темному с использованием горизонтального сглаживания (`-s`) и проверки минимальной длины темного участка (`-k`).
    * Найденные для нескольких последних строк края усредняются с помощью **медианного фильтра** (`-M`) для получения стабильного положения края для текущей строки.
* **Обрезка и сохранение:** Выполняется вертикальная обрезка сегмента. Затем для каждой строки данные копируются в выходной файл, начиная от **выровненного медианой левого края** для этой строки. Ширина выходного файла задается опцией `-w` или рассчитывается автоматически. Это позволяет компенсировать небольшие изгибы или наклон дорожки.

## Основные опции командной строки

```
Usage: ./split -i <in.pgm> -m <out_dir> -t <thresh> [-e] [options]

Required:
  -i <input.pgm>    : Input PGM file (P5 format, maxval <= 255)
  -m <pgm_dir>      : Output directory for segments
  -t <threshold>    : Absolute brightness threshold [0-255]

Mode Selection:
  -e                : Use V2 dynamic edge tracking logic (default: V1 static edge)

Common Optional:
  -w <width>        : Segment processing width (>=0, 0=auto ~avg*0.95, default: 0)
  -b                : Enable debug messages

V1 Logic Options (used if -e is NOT set):
  -s <smooth_V1>    : Initial column average smoothing window (>=1, default: 1)

V2 Logic Options (used if -e IS set):
  -S <smooth_I>     : Initial column average smoothing window (>=1, default: 1)
  -s <smooth_R>     : Row horizontal smoothing window for edge tracking (>=1, default: 1)
  -M <median_W>     : Median filter window size for edge tracking (odd >=3, default: 3)
  -d <delta>        : Edge search window delta (+/-, optional, default: final_width/2)
  -k <runlen>       : Min dark pixels run length to confirm edge (>=1, default: 1)

Примечание: В режиме V1 опции -S, -M, -d, -k игнорируются. Опция -s используется для начального сглаживания.
Примечание: В режиме V2 опция -S используется для начального сглаживания, а -s - для сглаживания строк.
```

## Разработка
Программа была написана с помощью ИИ `Gemini` на основе описанной мной логики и алгоритмов. Было решено отказаться от сложных алгоритмов компьютерного зрения (детального поиска краев, линий, фигур) в пользу более простого и прямолинейного подхода, достаточного для данной задачи.
