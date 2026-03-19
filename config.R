# ============================================================================
# 配置文件：NHANES HCF-4α 研究项目
# ============================================================================

# 自动检测项目根目录
if (!exists("PROJECT_ROOT")) {
  tryCatch({
    PROJECT_ROOT <- dirname(normalizePath(sys.frame(1)$ofile))
  }, error = function(e) {
    PROJECT_ROOT <- getwd()
  })
}

# ============================================================================
# 数据目录
# ============================================================================
RAW_DATA_DIR <- file.path(PROJECT_ROOT, "data", "raw")

# ============================================================================
# 处理后的数据目录
# ============================================================================
PROCESSED_DATA_DIR <- file.path(PROJECT_ROOT, "data", "processed")

L_DATA_DIR <- file.path(PROCESSED_DATA_DIR, "L_cycle")
P_DATA_DIR <- file.path(PROCESSED_DATA_DIR, "P_cycle")

# ============================================================================
# 输出目录（全部放在 processed 下）
# ============================================================================
L_RESULTS_DIR <- file.path(L_DATA_DIR, "results")
P_RESULTS_DIR <- file.path(P_DATA_DIR, "results")

# 论文特定结果目录（L周期）
PAPER1_RESULTS_DIR <- file.path(L_RESULTS_DIR, "paper1")
PAPER2_RESULTS_DIR <- file.path(L_RESULTS_DIR, "paper2")
PAPER3_RESULTS_DIR <- file.path(L_RESULTS_DIR, "paper3")
CROSS_CYCLE_RESULTS_DIR <- file.path(L_RESULTS_DIR, "cross_cycle")
LONGITUDINAL_RESULTS_DIR <- file.path(L_RESULTS_DIR, "longitudinal")
SUPPLEMENT_RESULTS_DIR <- file.path(L_RESULTS_DIR, "supplement")
SENSITIVITY_RESULTS_DIR <- file.path(L_RESULTS_DIR, "sensitivity")
SERIOUS_RESULTS_DIR <- file.path(L_RESULTS_DIR, "serious_research")

# 论文特定结果目录（P周期）
PAPER1_RESULTS_P_DIR <- file.path(P_RESULTS_DIR, "paper1")
PAPER2_RESULTS_P_DIR <- file.path(P_RESULTS_DIR, "paper2")
PAPER3_RESULTS_P_DIR <- file.path(P_RESULTS_DIR, "paper3")

# ============================================================================
# 日志目录（直接在根目录下）
# ============================================================================
LOGS_DIR <- file.path(PROJECT_ROOT, "logs")

# 创建所有目录
dirs <- c(RAW_DATA_DIR, PROCESSED_DATA_DIR, L_DATA_DIR, P_DATA_DIR,
          L_RESULTS_DIR, P_RESULTS_DIR,
          PAPER1_RESULTS_DIR, PAPER2_RESULTS_DIR, PAPER3_RESULTS_DIR,
          CROSS_CYCLE_RESULTS_DIR, LONGITUDINAL_RESULTS_DIR,
          SUPPLEMENT_RESULTS_DIR, SENSITIVITY_RESULTS_DIR, SERIOUS_RESULTS_DIR,
          PAPER1_RESULTS_P_DIR, PAPER2_RESULTS_P_DIR, PAPER3_RESULTS_P_DIR,
          LOGS_DIR)

for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat("✅ 配置加载成功\n")
cat("项目根目录:", PROJECT_ROOT, "\n")