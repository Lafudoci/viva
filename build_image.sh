#!/bin/bash
# ============================================================
# VIVA Docker Image 自動建置腳本
# 自動抓取最新 Git Tag 作為版號，建置 Docker Image 並匯出
# ============================================================

set -euo pipefail

# --- 設定 ---
IMAGE_NAME="viva"
OUTPUT_DIR="$HOME/docker-images"

# 切換到腳本所在目錄（專案根目錄）
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# --- 取得最新 Git Tag ---
LATEST_TAG=$(git describe --tags --abbrev=0 2>/dev/null)
if [ -z "$LATEST_TAG" ]; then
    echo "❌ 錯誤：找不到任何 Git Tag，請先建立 Tag 後再執行。"
    exit 1
fi

FULL_IMAGE="${IMAGE_NAME}:${LATEST_TAG}"
OUTPUT_FILE="${OUTPUT_DIR}/${IMAGE_NAME}-${LATEST_TAG}.tar"

echo "========================================"
echo "  VIVA Docker Image 自動建置"
echo "========================================"
echo "  專案目錄：${SCRIPT_DIR}"
echo "  最新 Tag：${LATEST_TAG}"
echo "  Image 名稱：${FULL_IMAGE}"
echo "  匯出路徑：${OUTPUT_FILE}"
echo "========================================"

# --- 建立輸出目錄 ---
mkdir -p "$OUTPUT_DIR"

# --- Docker Build ---
echo ""
echo "🔨 開始建置 Docker Image: ${FULL_IMAGE} ..."
echo ""
sudo docker build -t "$FULL_IMAGE" .

echo ""
echo "✅ Docker Image 建置完成：${FULL_IMAGE}"

# --- 匯出 Image ---
echo ""
echo "💾 匯出 Image 至 ${OUTPUT_FILE} ..."
sudo docker save -o "$OUTPUT_FILE" "$FULL_IMAGE"
sudo chown "$(id -u):$(id -g)" "$OUTPUT_FILE"

echo ""
echo "========================================"
echo "✅ 全部完成！"
echo "  Image：${FULL_IMAGE}"
echo "  檔案：${OUTPUT_FILE}"
echo "  大小：$(du -h "$OUTPUT_FILE" | cut -f1)"
echo "========================================"
