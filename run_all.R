#!/usr/bin/env Rscript

# ============================================================================
# 自动检测项目根目录
# ============================================================================
if (!exists("PROJECT_ROOT")) {
  tryCatch({
    PROJECT_ROOT <- dirname(normalizePath(sys.frame(1)$ofile))
  }, error = function(e) {
    PROJECT_ROOT <- getwd()
  })
}
setwd(PROJECT_ROOT)

source("config.R")

cat("\n", rep("=", 60), "\n", sep="")
cat("NHANES HCF-4α 研究 - 完整分析流程\n")
cat("开始时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 60), "\n\n", sep="")

start_time <- Sys.time()

# ============================================================================
# 第一阶段：数据清洗和基础变量构建（01-04）
# ============================================================================
cat("\n=== 第一阶段：数据清洗和基础变量构建 ===\n")

phase1_scripts <- c(
  "01_data_cleaning_L.R",
  "01_data_cleaning_P.R",
  "02_data_merge_L.R",
  "02_data_merge_P.R",
  "03_variable_basic_L.R",
  "03_variable_basic_P.R",
  "04_variable_biomarker_L.R",
  "04_variable_biomarker_P.R"
)

for (script in phase1_scripts) {
  if (file.exists(script)) {
    cat("\n--- 运行:", script, "---\n")
    source(script)
  } else {
    cat("⚠️ 文件不存在:", script, "\n")
  }
}

# ============================================================================
# 第二阶段：核心分析（08-22）- 依赖05-07生成的RDS文件
# ============================================================================
cat("\n=== 第二阶段：核心分析 ===\n")
cat("注：本阶段依赖05-07生成的RDS文件，这些文件已放置在 data/processed/ 目录下\n")

phase2_scripts <- c(
  "08_final_merge_L.R",
  "08_final_merge_P.R",
  "09_paper1_HCF_analysis_L.R",
  "09_paper1_HCF_analysis_P.R",
  "10_paper1_supplement_L.R",
  "10_paper1_supplement_P.R",
  "11_paper1_deep_analysis_L.R",
  "11_paper1_deep_analysis_P.R",
  "12_paper2_pathway_analysis_L.R",
  "12_paper2_pathway_analysis_P.R",
  "13_paper2_supplement_L.R",
  "13_paper2_supplement_P.R",
  "14_paper2_deep_analysis_L.R",
  "14_paper2_deep_analysis_P.R",
  "15_paper2_deep_analysis_continued_L.R",
  "15_paper2_deep_analysis_continued_P.R",
  "16_paper3_alpha_analysis_L.R",
  "16_paper3_alpha_analysis_P.R",
  "17_cross_cycle_validation_PL.R",
  "18_longitudinal_analysis_PL.R",
  "19_supplement_analysis_PL.R",
  "20_sensitivity_analysis_PL.R",
  "21_serious research.R",
  "22补充 等待审稿人要求.R"
)

for (script in phase2_scripts) {
  if (file.exists(script)) {
    cat("\n--- 运行:", script, "---\n")
    source(script)
  } else {
    cat("⚠️ 文件不存在:", script, "\n")
  }
}

# ============================================================================
# 完成
# ============================================================================
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n", rep("=", 60), "\n", sep="")
cat("✅ 所有分析完成！\n")
cat("完成时间:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("总耗时:", round(duration, 2), "分钟\n")
cat(rep("=", 60), "\n", sep="")