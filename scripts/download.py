import os

# 读取存好的链接文件
with open('../urls.txt', 'r') as f:
    urls = [line.strip() for line in f if line.strip()]

# 创建保存数据的目录
os.makedirs('../raw_data', exist_ok=True)

for url in urls:
    filename = url.split('/')[-1]
    print(f"正在下载: {filename}...")
    # 注意：Windows 需安装 wget，或在 Linux 环境运行
    os.system(f"wget -P ../raw_data/ {url}")
