#!/bin/bash
# ============================================================
# VIVA Docker Image 自動建置腳本
# 根據當前 Git 狀態（Tag 或 Commit）建置 Docker Image 並匯出
# ============================================================

set -euo pipefail

# --- 設定 ---
IMAGE_NAME="viva"
OUTPUT_DIR="$HOME/docker-images"

# 切換到腳本所在目錄（專案根目錄）
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# --- 取得當前 Git 版本描述 ---
# 若當前 commit 正好是 tag，則為 tag 名稱
# 若不是，則會呈現 v1.19.1-1-g487d2d9 這種格式
VERSION=$(git describe --tags --always 2>/dev/null || echo "latest")

FULL_IMAGE="${IMAGE_NAME}:${VERSION}"
OUTPUT_FILE="${OUTPUT_DIR}/${IMAGE_NAME}-${VERSION}.tar"

echo "========================================"
echo "  VIVA Docker Image 自動建置"
echo "========================================"
echo "  專案目錄：${SCRIPT_DIR}"
echo "  建置版本：${VERSION}"
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

# 複製 Wrapper Script
WRAPPER_FILE="${OUTPUT_DIR}/viva"
echo "📄 複製 Wrapper Script 至 ${WRAPPER_FILE} ..."
cp "$SCRIPT_DIR/viva" "$WRAPPER_FILE"
chmod +x "$WRAPPER_FILE"

echo ""
echo "========================================"
echo "✅ 全部完成！"
echo "  Image：${FULL_IMAGE}"
echo "  檔案：${OUTPUT_FILE}"
echo "  大小：$(du -h "$OUTPUT_FILE" | cut -f1)"
echo "========================================"
