# run_export.R

# 1. 加载依赖
library(hdf5r)
library(Seurat)
library(Matrix)

# 2. 加载模块
source("utils_write_matrix.R")
source("utils_write_metadata.R")
source("utils_write_components.R")
source("export_h5ad.R")

# 3. 假设你内存里已经有了 seurat_obj (刚刚转出来的那个)
if (!exists("seurat_obj")) stop("请先加载一个 seurat_obj")

# 4. 执行导出
output_file <- "export_test.h5ad"
seurat_to_h5ad(seurat_obj, output_file)

# 5. 自检 (可选：尝试用你的读取脚本再读回来，看是否一致)
# 这是一个完美的闭环测试！
message("\n>>> [Loop Test] 正在尝试读取刚刚导出的文件...")
source("run_conversion.R") # 确保这个脚本里调用的是 output_file
# 你可以临时修改 run_conversion.R 里的 target_file <- "export_test.h5ad"