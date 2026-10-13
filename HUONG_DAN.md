# Chạy pipeline phân tích 16S bằng Python và R

`run_pipeline.py` mở cửa sổ chọn dữ liệu. Bản Windows phát hành dùng PyInstaller
để đóng gói Python và giao diện; ứng dụng điều phối một tiến trình R;
phân tích vẫn do DADA2, phyloseq, vegan, DECIPHER, phangorn và ANCOMBC thực hiện.
Toàn bộ 17 script đã được sao chép vào thư mục `pipeline` cạnh file Python.

## 1. Chuẩn bị

- **Người dùng bản phát hành:** cài một lần bằng `Databio16S-Setup.exe` rồi mở
  lối tắt Databio 16S. Không cần cài Python, RStudio, Cutadapt thủ công hay sửa
  PATH. Lần mở đầu, trình thiết lập dò R, xác minh Cutadapt và cài package R
  còn thiếu. Cần Internet để tải package R.
- **R được cài riêng:** R/Bioconductor là môi trường khoa học lớn, độc lập và
  không được đóng gói trong ứng dụng. Nếu chưa có R, trình thiết lập hướng dẫn
  tải bản R chính thức; sau khi cài, chọn `Rscript.exe` nếu tự dò không tìm
  thấy. Ứng dụng chỉ dùng `Rscript.exe`; không yêu cầu RStudio.
- **Cutadapt:** được đóng gói trong bản Windows tạo bằng script phát hành. Khi
  chạy mã nguồn trực tiếp, cần chọn executable Cutadapt đã cài.
- **Nhà phát triển chạy mã nguồn:** cần Python 3.9 trở lên có Tcl/Tk. Giao diện
  không cần package Python bổ sung.
- FASTQ phải là **paired-end**, tên như `Sample01_R1.fastq.gz` và
  `Sample01_R2.fastq.gz`, hoặc `_R1_001/_R2_001`, `_1/_2`.
  Phần tên trước hậu tố là SampleID đầy đủ; mỗi SampleID chỉ có một cặp.
- SILVA phải là **training FASTA định dạng DADA2**, chẳng hạn
  `silva_nr99_v138.2_toGenus_trainset.fa.gz`, không phải FASTA bất kỳ.

### Tạo bộ cài Windows

Xem [packaging/README.md](./packaging/README.md). Trên máy đóng gói Windows
cần Python 3.12 và Inno Setup; chạy `.\packaging\build_windows.ps1` để tạo ZIP
di động và, khi tìm thấy Inno Setup, `dist\Databio16S-Setup.exe`. Mã nguồn hiện
không chứa bộ cài hay binary Cutadapt đã đóng gói; script tạo chúng từ các
dependency có phiên bản ghim.

## 2. Mở ứng dụng

Người dùng bản phát hành mở lối tắt **Databio 16S** sau khi cài Setup.exe.
Để chạy mã nguồn trong lúc phát triển, dùng `start_pipeline.cmd` hoặc terminal:

```powershell
cd "C:\Users\Administrator\Downloads\databio_16S_python_R"
py -3 run_pipeline.py
```

1. Ở lần chạy đầu, trình thiết lập kiểm tra Rscript và Cutadapt rồi cài package
   R còn thiếu. Chọn ngôn ngữ nếu cần; sau khi hoàn tất ứng dụng mở giao diện.
2. Chọn **Ngôn ngữ / Language** ở góc trên bên phải để dùng tiếng Việt hoặc
   tiếng Anh. Lựa chọn này được lưu trong `%LOCALAPPDATA%\Databio16S`.
3. Bấm **Chọn file…** và dùng Ctrl/Shift để chọn tất cả R1/R2; hoặc bấm
   **Chọn thư mục…** để tìm FASTQ trong thư mục và thư mục con.
   Các thư mục kết quả có tên `result`, `results`, `runs`, `filtered`,
   `primer_trimmed` và thư mục đầu ra đã chọn sẽ được bỏ qua khi quét.
   Nếu dữ liệu gốc có tên thư mục như vậy, dùng nút chọn file trực tiếp.
4. Chọn **metadata CSV**, **SILVA trainset**, và **thư mục lưu kết quả**.
5. Chọn cột nhóm và nhóm tham chiếu. Primer, ASV min/max và ngưỡng Q được
   lấy mặc định từ cấu hình đã gửi, cần kiểm tra theo thí nghiệm của bạn.
   Các thông số khác được giữ trong `pipeline/00_parameters.R`.
6. Bấm **Kiểm tra môi trường / dữ liệu**: kiểm tra package R, Cutadapt,
   ghép cặp FASTQ và metadata; chưa xử lý reads hay kiểm chứng reference.
7. Khi kiểm tra thành công, bấm **Chạy phân tích**. Theo dõi log bên dưới.
   Giao diện giữ phản hồi trong khi R chạy. Có thể bấm **Dừng** để ngừng
   cả tiến trình R và Cutadapt con; kết quả bị dừng được đánh dấu chưa hoàn chỉnh.

Ứng dụng lưu thiết lập trong `%LOCALAPPDATA%\Databio16S\settings.json`, do đó
vẫn lưu được khi cài trong Program Files. Khi chạy mã nguồn, giữ
`run_pipeline.py`, `app_setup.py`, `run_pipeline.R`, `install_packages.R` và
thư mục `pipeline` cùng nhau.

## 3. Metadata mẫu

CSV UTF-8, phân cách bằng dấu phẩy, cột đầu là SampleID:

```csv
SampleID,Group
Sample01,Control
Sample02,Control
Sample03,Treatment
Sample04,Treatment
```

SampleID phải khớp tên FASTQ, không bao gồm `_R1`/`_R2` hay phần mở rộng.
Được phép dùng SampleID có số 0 ở đầu như `001`; số 0 được giữ nguyên.
Tên cột nhóm dùng chữ, số, dấu gạch dưới (ví dụ `Group`, `Treatment`).
Pipeline này thực hiện PERMANOVA và kiểm tra phân tán, nên giao diện yêu cầu
ít nhất hai nhóm và hai mẫu mỗi nhóm. Đây là điều kiện chạy kỹ thuật tối thiểu,
không thay thế thiết kế thí nghiệm hay tính toán số mẫu.
ANCOM-BC có thể tắt bằng checkbox; các phép so sánh beta vẫn chạy.

## 4. Kết quả

Mỗi lần chạy tạo một thư mục riêng, không ghi đè lần chạy trước:

```text
Documents/Databio16S/runs/analysis_YYYYMMDD_HHMMSS_microseconds/
  config.json                 # Đường dẫn và thông số đã chọn
  input_manifest.csv          # SampleID và file R1/R2
  pipeline.log                # Toàn bộ stdout/stderr của R
  status.json                 # Mã thoát và trạng thái dừng
  result/
    figures/                  # PDF/PNG
    tables/                   # CSV và FASTA
    rds/                      # Đối tượng phân tích R
    logs/                     # Log Cutadapt, sessionInfo.txt, ngưỡng trim
    primer_trimmed/
    filtered/
```

Nút **Mở kết quả** mở thư mục lần chạy vừa xong. Lần kiểm tra tạo thư mục
`check_...`, không tạo kết quả phân tích. Mỗi lần bấm Chạy là một lần chạy mới;
giao diện chưa hỗ trợ tiếp tục một thư mục lần chạy cũ.
Không xóa FASTQ đầu vào. Chọn file gốc nếu muốn pipeline tự cắt primer.
Nếu bạn đã cắt primer của dữ liệu bằng Cutadapt trước đó, bước 03 vẫn cắt lại
và có `--discard-untrimmed`; bộ mã hiện tại chưa có chế độ bỏ qua cắt primer.

## 5. Thay đổi so với các script đã gửi

- Bản gốc ngoài thư mục ứng dụng được giữ nguyên.
- Python ghi cấu hình JSON, R nạp bước 00 rồi ghi đè các lựa chọn từ giao diện
  trước khi chạy bước 01 tới 16 trong cùng một môi trường R.
- Bước 02 chỉ lấy FASTQ đã chọn, tránh quét nhầm output vào đầu vào.
- Metadata giữ SampleID dưới dạng ký tự và hỗ trợ UTF-8 BOM.
- Bước 03 gọi executable Cutadapt đã chọn, quote từng đối số để xử lý đường
  dẫn có dấu cách, và dừng ngay khi Cutadapt thất bại hoặc thiếu output.
- File cài đặt liệt kê `phangorn` trong CRAN, bổ sung `jsonlite` và `rlang`;
  không chứa đường dẫn Cutadapt cố định.
- Kiểm tra trước các package, nhóm, reference group và file đầu vào.
- Phương pháp thống kê và lọc sinh học còn lại theo script gốc. Ví dụ danh sách
  `contaminant_genera` được giữ nguyên; cần rà soát theo nghiên cứu của bạn.
  Bước chọn truncation đọc toàn bộ FASTQ vào RAM theo script gốc; dữ liệu lớn
  có thể cần máy nhiều RAM. Số luồng mặc định là 1.

## 6. Khi gặp lỗi

- **Không tìm thấy Rscript**: cài R từ trang chính thức, sau đó dùng nút
  **Thiết lập thành phần** để chọn `Rscript.exe`. Không cần RStudio hay sửa PATH.
- **Thiếu package R / cài package thất bại**: mở **Thiết lập thành phần** để
  thử lại; kiểm tra Internet và dung lượng ổ đĩa. Package được cài vào thư viện
  người dùng trong R.
- **Không tìm thấy Cutadapt trong bản phát hành**: cài lại bản chính thức hoặc
  dùng **Thiết lập thành phần** để chọn executable hợp lệ. Bản chạy mã nguồn
  cần Cutadapt được cài riêng.
- **SampleID không khớp**: đối chiếu `input_manifest.csv` với cột đầu metadata.
- **Không ghép được R1/R2**: kiểm tra tên file, không chọn trùng mẫu giữa nhiều
  lane hoặc nhiều thư mục; đặt các cặp theo cùng cấu trúc thư mục.
- **Cutadapt lỗi**: xem `pipeline.log` và `result/logs/cutadapt_primer_removal.log`.
- **Mất nhiều reads / quá ít ASV**: kiểm tra primer, chất lượng, ngưỡng trim,
  chiều dài ASV, đủ overlap để merge và database phân loại; không chỉ tăng
  ngưỡng một cách tùy ý để pipeline chạy qua.
- Không mở hai cửa sổ ứng dụng để phân tích cùng lúc nếu máy không đủ RAM.

## 7. Kiểm thử mã

```powershell
py -3 -m unittest discover -s tests -v
```

Các test kiểm tra dò dependency, chọn Cutadapt, ghép cặp FASTQ, metadata,
thông báo lỗi và lưu trạng thái tiến trình. Chúng không chứng minh tính đúng
đắn khoa học trên dữ liệu thực hoặc thay thế việc kiểm thử bộ cài Windows.
Xem [VALIDATION.md](./VALIDATION.md) để biết phần nào đã và chưa được kiểm tra.
