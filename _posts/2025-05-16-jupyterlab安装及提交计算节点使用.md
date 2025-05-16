---
layout: post
title: "jupyterlab安装及提交计算节点使用"
date: 2025-05-16
description: "如何安装jupyter-lab，提交计算节点使用，如何在jupyter中使用R"

tag: 生信经验
---   



## 使用

conda activate scpEnv

sbatch run_jupyterLab_slum.sh，获取返回的日志文件中的计算节点和端口号，通过ssh隧道链接计算节点使用。

添加R kernel，在R中运行 IRkernel::installspec(name = 'R4.0.4', displayname = 'R4.0.4')

直接可以在base中使用

## 使用问题

可以直接在Console中安装包，如果用conda安装的包无法载入，可以通过找到配置文件，修改里面的env 对应的路径。

## **1.安装**

在服务器上安装 Jupyter Lab

```bash
# install
conda install -c conda-forge jupyterlab
# start
jupyter-lab
# 常带参数
jupyter-lab --no-browser --port 8889  # 不自动打开浏览器 + 指定 8889 端口
# 服务器中的 Jupyter Lab 需要常驻，选用了终端复用工具 tmux。
# 远程访问
# 生成配置文件
jupyter-lab --generate-config
jupyter server password
c.ServerApp.root_dir = '/icto/user/yc47650/data-home/analysisProject'
```

------

## **2.远程访问-修改config**

~/.jupyter/jupyter_lab_config.py 中的以下内容：

```bash
c.ServerApp.allow_remote_access = True
c.ServerApp.ip = '*'

# 根据个人需要修改端口号（默认 8888）
c.ServerApp.port = 8889
```

## 3. 端口转发

ssh隧道，跳板机，端口号保持一致

![image.png](https://prod-files-secure.s3.us-west-2.amazonaws.com/7b662f88-2133-4bff-abef-45290bc9e34a/28a35c19-bc86-4307-9d2e-33c67af6f84d/image.png)

## 4. 本地启动

此时，在本地电脑的浏览器地址栏中输入 `localhost:8887`，即可打开 Jupyter Lab。然后输入taken，在服务器启动时显示。

## 5. 安装内核

```
# 安装 ipykernel
$ conda install ipykernel

# 将环境写入notebook的kernel中
$ python -m ipykernel install --user --name 环境名称 --display-name "在jupyter中显示的环境名称"
# 比如已经有个叫 `gpu` 的 conda 环境，希望在 Jupyter 的内核选择中显示为 `deeplearning`
$ python -m ipykernel install --user --name gpu --display-name "deeplearning"

# 如果不带 --display-name 参数，则与与环境名同名，在上述的例子里就是 `gpu`
```

**安装内核bug**

```bash
# 已经正确安装了内核，但是一直无法载入包
jupyter kernelspec list # 查看内核路径
# 进入路径直接修改json文件，可以改环境和对应的内核名
 kernel.json 
```

## 6. R-kernel jupyter运行R

在R中安装IRkernel

conda install conda-forge::r-irkernel

然后在R中运行IRkernel::installspec(name = 'scpEnv', displayname = 'scpR')让jupyter可以找到R

## 连接远程服务器

https://jinsblog.com/posts/jnb-remote-server/

这里的YourHost就是ssh连接的HPCC，主机名，

节点sbatch以下代码，然后把节点服务器地址拷贝过来，输入对应端口号就可以了。在本地网页输入http://127.0.0.1:8899/lab



```shell
#!/bin/bash

### HIC 2024

#SBATCH -J jupyter
#SBATCH -n 24                              # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
#SBATCH -t 1-00:00:0                         # Runtime in D-HH:MM format
#SBATCH -p fhs-fast                           # Partition to run in
#SBATCH --mem=100G                        # Memory total in MB (for all cores)
#SBATCH --output=bc1_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH --error=bc1_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL,END                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=yc47650@um.edu.mo

export TMPDIR=/scratch2/$USER/$SLURM_JOB_ID
mkdir -p $TMPDIR

# get tunneling info
XDG_RUNTIME_DIR=""
port=$(shuf -i8000-9999 -n1)
node=$(hostname -s)
user=$(whoami)
cluster=$(hostname -f | awk -F"." '{print $2}')

# print tunneling instructions jupyter-log
echo -e "
   MacOS or linux terminal command to create your ssh tunnel
   ssh YourHost -N -L 8890:${node}:${port}

   Use a Browser on your local machine to go to:
   localhost:${port}  (prefix w/ https:// if using password)
   "

# load modules or conda environments here

source ~/software/miniconda3/etc/profile.d/conda.sh
conda activate scpEnv

cd /icto/user/yc47650/data-home/analysisProject # where you run this job
jupyter lab --no-browser --port=${port} --ip=${node}
# jupyter-notebook --no-browser --port=${port} --ip=${node}
```



