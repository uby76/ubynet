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

### 3. 基于强度变化的反应匹配

```r
match_res <- match_reactions_by_intensity(
  file1 = "inflow.csv",
  file2 = "outflow.csv",
  reaction_delta_file = "reaction_deltas.csv",
  out_dir = "results",
  intensity_ratio_lower = 0.5,
  intensity_ratio_upper = 2,
  mass_tolerance = 0.005
)
```

### 4. 基于质量差的反应匹配

```r
match_res2 <- match_reactions_by_mass_difference(
  file1 = "inflow.csv",
  file2 = "outflow.csv",
  reaction_delta_file = "reaction_deltas.csv",
  out_dir = "results",
  mass_tolerance = 0.005
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

* **功能**：基于两个样本的强度变化 + 质量差，匹配潜在反应。
* **输入**：两个 CSV 文件 + 反应定义表
* **输出**：反应网络边表、反应统计摘要

---

### `match_reactions_by_mass_difference()`

* **功能**：基于质量差匹配反应（不考虑强度比）。
* **输入**：两个 CSV 文件 + 反应定义表
* **输出**：反应网络边表、反应统计摘要

---

### `build_mass_pmd_network()`

* **功能**：基于前体与产物集合构建 PMD 网络。
* **输入**：前体数据表、产物数据表、反应定义表
* **输出**：网络边表（包含 Source, Target, Reaction, Mass Error 等）

---

### `transformation_analysis()`

* **功能**：对 PMD 网络进行进一步分析（如路径、特征统计）。
* **输入**：网络边表
* **输出**：分析结果

---

### `classification_functions()`

* **功能**：对分子/反应进行分类与标签。
* **输入**：分子或反应数据
* **输出**：分类结果

---

### `compare_molforms()`

* **功能**：比较两个分子式的元素差异。
* **输入**：两个分子式（字符串）
* **输出**：元素组成差异

---

## 📊 函数总览表

| 函数                                   | 功能         | 输入            | 输出                     |
| ------------------------------------ | ---------- | ------------- | ---------------------- |
| `merge_mass_intensity`               | 合并质量-强度数据  | CSV 文件夹       | 合并 mass-intensity 表    |
| `merge_molform_intensity`            | 合并分子式-强度数据 | CSV 文件夹       | 合并 molform-intensity 表 |
| `match_reactions_by_intensity`       | 基于强度变化匹配反应 | 2 个 CSV + 反应表 | 边表 + 反应摘要              |
| `match_reactions_by_mass_difference` | 基于质量差匹配反应  | 2 个 CSV + 反应表 | 边表 + 反应摘要              |
| `build_mass_pmd_network`             | 构建 PMD 网络  | 前体 + 产物 + 反应表 | 网络边表                   |
| `transformation_analysis`            | 分子转化网络分析   | 网络边表          | 分析结果                   |
| `classification_functions`           | 功能分类       | 分子或反应数据       | 分类结果                   |
| `compare_molforms`                   | 分子式差异比较    | 两个分子式         | 元素差异                   |

