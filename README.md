# ubynet

**ubynet** 是一个用于 **分子数据分析与转化网络构建** 的 R 包。  
它提供了从分子数据合并、质量差反应匹配，到构建基于成对质量差 (PMD) 的分子转化网络的完整工具链。  
适合处理 FTICR-MS 等超高分辨质谱数据的反应分析与分子转化研究。

---

## 🔧 安装

确保已安装 [`remotes`](https://cran.r-project.org/package=remotes)：

```r
install.packages("remotes")
````

从 GitHub 安装 `ubynet`：

```r
remotes::install_github("uby76/ubynet")
```

---

## 🚀 使用示例

### 1. 合并质量-强度数据

```r
res <- merge_mass_intensity(
  dir_path = "data/csv_files",
  output_mass_intensity = "out/mass_int.csv",
  output_mass_elements = "out/mass_el.csv",
  output_meta_file = "out/meta.csv"
)
```

### 2. 合并分子式-强度数据

```r
res2 <- merge_molform_intensity(
  dir_path = "data/csv_files",
  output_molform_intensity = "out/molint.csv",
  output_molform_elements = "out/molel.csv"
)
```

### 3. 基于分子式变化的反应匹配

```r
match_res <- match_reactions_by_intensity(
  file1 = "inflow.csv",
  file2 = "outflow.csv",
  reaction_delta_file = "reaction_deltas.csv",
  out_dir = "results"
)
```

### 4. 基于质量差的反应匹配

```r
match_res2 <- match_reactions_by_mass_difference(
  file1 = "inflow.csv",
  file2 = "outflow.csv",
  reaction_delta_file = "reaction_deltas.csv",
  out_dir = "results",
  mass_tolerance = 0.005 #设置容差
)
```


---

### 5. 基于已知的分子式和MASS计算的PMD反应网络
## 参考文献：https://www.nature.com/articles/s41467-020-19989-y

```r
edges <- build_mass_pmd_network(
  mol_file = "MS_MolInfor1.csv",
  trans_file = "Transformation_Database_07-2020.csv",
  error_term = 0.00001,
  output_dir = "MS_MolInfor2"
)

```
---

## 📖 函数说明

### `merge_mass_intensity()`

* **功能**：合并质量与强度数据，生成统一的质量-强度表。
* **输入**：CSV 文件目录
* **输出**：合并后的 mass-intensity 数据表 + 元数据

---

### `merge_molform_intensity()`

* **功能**：合并分子式与强度数据。
* **输入**：CSV 文件目录
* **输出**：合并后的 molform-intensity 数据表 + 元数据

---

### `match_reactions_by_intensity()`

* **功能**：基于两个样本的强度变化 + 分子式变化，匹配潜在反应。
* **输入**：两个 CSV 文件 + 反应定义表
* **输出**：反应网络边表、反应统计摘要

---

### `match_reactions_by_mass_difference()`

* **功能**：基于质量差匹配反应考虑+两个样本的强度变化。
* **输入**：两个 CSV 文件 + 反应定义表
* **输出**：反应网络边表、反应统计摘要

---

### `build_mass_pmd_network()`
* **功能**：基于分子信息文件和反应数据库，构建成对质量差 (PMD) 分子转化网络。  
* **输入**：  
  - `mol_file`：分子信息文件 (CSV)，需包含分子式和质量信息  
  - `trans_file`：反应数据库文件 (CSV)，需包含 `reaction` 和 `mass_difference` 列  
  - `error_term`：质量差匹配容差 (默认 1e-5 Da)  
  - `output_dir`：结果输出目录  
* **输出**：  
  - 在 `output_dir` 中生成 PMD 网络边表与相关结果文件  
  - 返回构建好的网络边表  
 


## 📊 函数总览表

| 函数                                   | 功能         | 输入            | 输出                     |
| ------------------------------------ | ---------- | ------------- | ---------------------- |
| `merge_mass_intensity`               | 合并质量-强度数据  | CSV 文件夹       | 合并 mass-intensity 表    |
| `merge_molform_intensity`            | 合并分子式-强度数据 | CSV 文件夹       | 合并 molform-intensity 表 |
| `match_reactions_by_intensity`       | 基于强度变化匹配反应 | 2 个 CSV + 反应表 | 边表 + 反应摘要              |
| `match_reactions_by_mass_difference` | 基于质量差匹配反应  | 2 个 CSV + 反应表 | 边表 + 反应摘要              |
| `build_mass_pmd_network()` | PMD网络  | 样本数据和已知的反应数据库 | 边表 + 反应摘要              |
