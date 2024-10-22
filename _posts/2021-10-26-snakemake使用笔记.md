layout: post
title: "snakemake使用笔记"
date: 2021-10- 26
description: "snakemake是一个基于python3的任务流程工具，执行流程脚本简单易懂，支持断点运行、并行运算、内存控制、cpu核心控制、流程控制等特点"

### snakemake命令与规则

#### 执行逻辑

snakemake并没有准备专门的顺序或依赖的语法关键词, 流程的串联完全依靠输入和输出文件的依赖关系自动完成. 由于snakemake默认会去解析并完成第一个`rule`, 因此官方文档推荐的使用方式是创建一个名称为”all”的`rule`, 然后在这个`rule`的`input`中(注意不是`output`)指定所有的最终文件. 然后程序会发现结果文件不存在, 便在文件中寻找到某条`rule`的`output`符合文件模式, 然后继续解析这个`rule`需要要什么作为`input`, 如果`input`(即执行依赖)已满足, 则会开始运行该规则下的命令/脚本, 如果不满足, 则会继续向前查找, 直到找到源头或者找不到源头而报错为止

#### rule中引用变量

```python
path="path/to/input"

rule samtools_flt:
    input: 
        samples = expand("{sample_id}.{fq}.fasta", sample_id=['sam1', 'sam2'], fq=['fq1', 'fq2'])
    output:
        "test.flt.bam"
    shell:
        "samtools sort {path}/{{input}} > {{output}}".format(path=path)
```

#### rule

​        `rule`是Snakefile中最主要的部分。每一个rule定义了pipeline流程中运行的一步，每一个rule都可以当作一个shell脚本来处理，一般主要包括 `input`、`output`、`shell` 3个部分。

1. `rule all`

   ​		不同于其他的rule，在`rule all`里面一般不会去定义要执行的命令，他一般用来定义最后的输出结果文件，rule all中的input是流程的最终目标文件。这个rule all 必须要有，不然整个pipeline不会运行。 

   ```python
   rule all:
       input:
            #"merged.txt"  取消注释后，则能正常输出文件
   rule concat:
       input:
           expand("{file}.txt", file=["hello", "world"])
       output:
           "merge.txt"
       shell:
           "cat {input} > {output}"
   ```

2. `wildcards`

   用来获取**通配符**匹配到的部分，例如对于通配符`"/file.{group}.txt"`匹配到文件`101/file.A.txt`，则`{wildcards.dataset}`就是101，`{wildcards.group}`就是A。

3. `threads`

   通过在rule里面指定`threads`参数来指定分配给程序的线程数，eg`threads: 8`。

4. `resources`

   可用来指定程序运行的内存，eg. `resources: mem_mb=800`。

5. `message`

   使用`message`参数可以指定每运行到一个rule时，在终端中给出提示信息，eg.`message: "starting mapping ..."`。

6. `priority`

   可用来指定程序运行的优先级，默认为0，eg.`priority: 20`。

7. `log` 

   用来指定生成的日志文件，eg.`log: "logs/concat.log"`。

8. `params`

   指定程序运行的参数，eg.`params: cat="-n"`,调用方法为`{params.cat}`。

9. `run`

   在`run`的缩进区域里面可以输入并执行python代码。

10. `scripts`

​		用来执行指定脚本。`script: 'scripts/script.py'; 'script/script.R'`

11. `temp`

    通过`temp`方法可以在所有`rule`运行完后删除指定的中间文件，如果下一步成功运行，那么就删掉这一步的输出文件，这样去掉没有必要保留的中间文件，就会为我们省掉一定的磁盘空间，大规模分析时候非常方便。 

    ```python
    fw = temp(WORKING_DIR + "fastq/{sample}_1.fastq"),
    rev= temp(WORKING_DIR + "fastq/{sample}_2.fastq") 
    ```

12. `protected`

    用来指定某些中间文件是需要保留的，eg.`output: protected("f1.bam")`。这里protected的意思就是输出文件很重要，我们想把访问权限改成只读，以防止误删和修改。 

    ```python
     bams  = protected(WORKING_DIR + "mApped/{sample}.bam") 
    ```

13. `ancient`

    重复运行执行某个Snakefile时，snakemake会通过比较输入文件的时间戳是否更改(比原来的新)来决定是否重新执行程序生成文件，使用ancient方法可以强制使得结果文件

14. Rule Dependencies

    可通过快捷方式指定前一个rule的输出文件为此rule的输入文件

    ```python
    rule a:
        input:  "path/to/input"
        output: "path/to/output"
        shell:  ...
    rule b:
        input:  rules.a.output   #直接通过rules.a.output 指定rule a的输出
        output: "path/to/output/of/b"
        shell:  ...
    ```

15. `report`

    使用snakemake定义的`report`函数可以方便的将结果嵌入到一个HTML文件中进行查看。

    ```python
    rule report:
        input:
            "calls/all.vcf"
        output:
            "report.html"
        run:
            from snakemake.utils import report
            with open(input[0]) as vcf:
                n_calls = sum(1 for l in vcf if not l.startswith("#"))
            report("""
            An example variant calling workflow
            ===================================
            Reads were mapped to the Yeast
            reference genome and variants were called jointly with
            SAMtools/BCFtools.
            This resulted in {n_calls} variants (see Table T1_).
            """, output[0], T1=input[0])
    ```

16. expand

    expand是一个特殊函数，它的目的是用通配符，将所有的SRA ID 传入。

    ```python
    qcfile = expand(RESULT_DIR + "fastp/{sample}.html",sample=SAMPLES) 
    ```

17. 自动创建文件夹

    在写Snakemake的主体代码部分时候我们只给出了文件夹链接，并不需要自己创建这些文件

#### 配置文件

​		每计算一次数据都要重写一次Snakefile有时可能会显得有些繁琐，我们可以将那些改动写入配置文件，使用相同流程计算时，将输入文件的文件名写入配置文件然后通过Snakefile读入即可。配置文件有两种书写格式——json和[yaml](https://link.jianshu.com?t=http%3A%2F%2Fwww.ruanyifeng.com%2Fblog%2F2016%2F07%2Fyaml.html)。在Snakefile中读入配置文件使用如下方式：

```python
configfile: "path/to/config.json" 
configfile: "path/to/config.yaml"
# 也可直接在执行snakemake命令时指定配置
$ snakemake --config yourparam=1.5
# config.yaml内容为：
samples:
    A: data/samples/A.fastq
    B: data/samples/B.fastq
# snakemake调用配置文件参数
configfile: "config.yaml"
...
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["smaples"])
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
# 虽然samples是一个字典，但是展开的时候，只会使用他们的key值部分。
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

#### snakemake执行

​		一般讲所有的参数配置写入Snakefile后直接在Snakefile所在路径执行`snakemake`命令即可开始执行流程任务。一些常用的参数：

```python
--snakefile, -s 指定Snakefile，否则是当前目录下的Snakefile
--dryrun, -n  不真正执行，一般用来查看Snakefile是否有错
--printshellcmds, -p   输出要执行的shell命令
--reason, -r  输出每条rule执行的原因,默认FALSE
--cores, --jobs, -j  指定运行的核数，若不指定，则使用最大的核数
--force, -f 重新运行第一条rule或指定的rule
--forceall, -F 重新运行所有的rule，不管是否已经有输出结果
--forcerun, -R 重新执行Snakefile，当更新了rule时候使用此命令
#一些可视化命令
$ snakemake --dag | dot -Tpdf > dag.pdf
#集群投递
snakemake --cluster "qsub -V -cwd -q 节点队列" -j 10
# --cluster /-c CMD: 集群运行指令
# qusb -V -cwd -q， 表示输出当前环境变量(-V),在当前目录下运行(-cwd), 投递到指定的队列(-q), 如果不指定则使用任何可用队列
# --local-cores N: 在每个集群中最多并行N核
# --cluster-config/-u FILE: 集群配置文件
```

#### snakemake运行逻辑

​		当Snakemake收到`rule_all`定义的需要输出文件后， 首先看这些输出文件是否存在，如果不存在是哪个`rule `产生了这些`output`，然后去运行相应的`rule`。


​		例如,我们在`rule all`里写了需要输出bam文件，那么`Snakemake `会先去看是`rule hisat_mApping` 产生了BAM，所以要先去运行`rule hisat_mApping`。在运行`rule hisat_mApping`的时候，它会发现`rule hisat_mApping`需要另外的input `trimed fastq`。

​		所以下一步是顺藤摸瓜，去看trimed fastq是哪个rule产生的，如此循环往复，直到回溯到第一步的SRA ID文件。运行的逻辑就像剥洋葱一样倒着把pipeline串起来。 如果是因为某一个原因pipeline中断了，假如运行到第二步 trimed fastq已经产生了，但是断电了。下次重新启动的时候，因为Snakemake是回溯后针对性的来运行，所以回溯到第二步trimed fastq，发现了这个文件，就不会在往前回溯了，保证了从断点处继续跑。

​		它执行分为三个阶段:

  1. 在初始化阶段，工作流程会被解析，所有规则都会被实例化

 2. 在DAG阶段，也就是生成有向无环图，确定依赖关系的时候，所有的通配名部分都会被真正的文件名代替。

3. 在调度阶段，DAG的任务按照顺序执行。

  ​		也就是说在初始化阶段，我们是无法获知通配符所指代的具体文件名，必须要等到第二阶段，才会有wildcards变量出现。也就是说之前的出错的原因都是因为第一个阶段没通过。这个时候就需要输入函数推迟文件名的确定，可以用Python的匿名函数，也可以是普通的函数

### snakemake规则

1. 确定运行步骤，主要通过输出文件的确定来匹配，这两有两种方式一种是直接运行snakemake时指定输出；另一种是通过rule all，如果在调用snakemake时没有指定目标文件，snakemake将尝试执行snakefile中的第一条规则，所以根据这个属性一般snakefile的第一个规则为rule all，rule all中的input规则仅仅是指示应该收集哪些结果。

   ```python
   snakemake --use-conda plots/quals.svg --cores 1
   ```

2. Report--HTML格式报告输出

   ```python
   snakemake --report report.html
   ```

   将输出运行实时统计，工作流程可视化拓扑图，使用的软件和数据信息，此外也可以标注任何工作流的输出文件到报告的结论中。例如将输出文件Output的'plots/quals.sng'写入报告的方法，只需要将命令替换为report("plots/quals.svg", caption = 'report/calling.rst')

3. snakefile文件是包含指示如何从输入文件出发创建输出文件。

#### 配置工作路径

1. snakefile中的所有路径都是相对于snakemake的目录进行解释的。可以通过在snakefile文件中指定一个workdir来覆盖此行为。

   ```
   workdir: "path/to/workdir"
   ```

2. 通常最好只通过命令行设置工作目录，因为上面的方法限制了流程的可移植性。

#### wildcard---将一个规则应用于大量数据集，自动化解析多个命名

1. 例子

   定义两个通配符wildcards，`dataset`和`group`，通配符用于取代正则表达`.+`。如果规则的输出与请求的文件匹配，通配符匹配的子字符串将传播到输入文件和变量通配符，这里也将在shell命令中使用这些通配符。

   ```python
   rule complex_conversion:
       input:
           "/inputfile"
       output:
           "/file.{group}.txt"
       shell:
           "somecommand --group {wildcards.group} < {input} > {output}"
   # shell中不能直接访问通配符对象里的内容，需要调用对象
   ```

   例如，对于识别的`101/file.A.txt`文件，snakemake将根据识别文件的结果设置`dataset=101`和`group=A`.

   **注意：**

   1. 在同一个规则中`input`和`output`中使用的通配符名称必须一致，但是不同规则中可以使用不同的命名来匹配同一个文件。

   2. 在一个文件中使用多个通配符可能会造成歧义

      如.{group}.txt`规则来匹配`101.B.norm.txt`,将不清楚是`dataset=101.B`还是`dataset=101`

   3. 对于2的问题，可以考虑给一个正则表达，我们可以严格限制`datast`通配符的组成

      ```python
      output: "{dataset,\d+}.{group}.txt" # 限制dataset为数值
      # 第二种方式，通过wild_constraints
      rule complex_conversion:
          input:
              "{dataset}/inputfile"
          output:
              "{dataset}/file.{group}.txt"
          wildcard_constraints:
              dataset="\d+"
          shell:
              "somecommand --group {wildcards.group}  < {input}  > {output}"
      # 第三种方式，定义全局通配符限制规则
      wildcard_constraints:
          dataset="\d+"
      
      rule a:
          ...
      ```

#### 聚合体

1. 输入文件可以是python的列表，允许简单的对参数或者样本聚合

   ```python
   rule aggregate:
       input:
           ["{dataset}/a.txt".format(dataset=dataset) for dataset in DATASETS]
       output:
           "aggregated.txt"
       shell:
           ...
   ```

2. expand函数

   ```python
   rule aggregate:
       input:
           expand("{dataset}/a.txt", dataset=DATASETS)
       output:
           "aggregated.txt"
       shell:
           ...
   ```

   **注意**，这里`dataset`不是一个通配符

   也可以在expand函数中使用通配符

   ```
   expand("{{dataset}/a.{ext}}", ext=FORMATS)
   ```

#### 系统资源控制

1. threads, 设置snakemake全局核心对象`workflow.cores`

2. 资源Resources

   ```python
   rule a:
       input:     ...
       output:    ...
       resources:
           mem_mb=100
       shell:
           "..."
   ```

#### priorities优先级

1. 任何规则的默认优先级是0，如果一个规则有更高的优先级，调度程序将优先于所有准备同时执行但至少没有相同优先级的规则。

#### 日志文件log-files

1. 日志文件是写入执行信息

2. 日志文件可以被其他规则作为输入使用，但是它与`output`区别是，它不会在程勋执行错误时被删除

   ```python
   rule abc:
       input: "input.txt"
       output: "output.txt"
       log: "logs/abc.log"
       shell: "somecommand --log {log} {input} {output}"
   ```

3. 对于没有显式`log`参数的程序可以使用`2>{log}`的方式重定向

#### rule中非文件的参数

1. 有时，希望定义与规则体分开的某些参数。可以使用`params`

   ```python
   rule:
       input:
           ...
       params:
           prefix="somedir/{sample}"
       output:
           "somedir/{sample}.csv"
       shell:
           "somecommand -o {params.prefix}"
   ```

#### flag files

1. 有时在没有文件依赖的情况下，需要强制一些规则的执行。这可以通过“touch”空文件来实现，这些空文件表示某个任务已完成。

   ```python
   rule all:
       input: "mytask.done"
   
   rule mytask:
       output: touch("mytask.done")
       shell: "mycommand ..."
   ```

#### 流程运行完成提示

1. 流程运行完成或错误提示

   ```
   onsuccess:
       print("Workflow finished, no error")
   
   onerror:
       print("An error occurred")
       shell("mail -s "an error occurred" youremail@provider.com < {log}")
   ```

#### rule dependencies

1. 在snakefile中其他规则的输出可以被引用

   ```python
   rule a:
       input:  "path/to/input"
       output: "path/to/output"
       shell:  ...
   
   rule b:
       input:  rules.a.output
       output: "path/to/output/of/b"
       shell:  ...
           
   rule a:
       input:  "path/to/input"
       output: a = "path/to/output", b = "path/to/output2"
       shell:  ...
   
   rule b:
       input:  rules.a.output.a
       output: "path/to/output/of/b"
       shell:  ...
   ```

#### 操作模糊规则

1. 当两个规则可以生成相同的输出文件时，如果没有额外的指导，snakemake无法决定使用哪一个。规则的顺序并不打算以正确的执行顺序引入规则(这完全取决于您使用的输入和输出文件的名称)，它只是帮助snakemake决定当多个规则可以创建相同的输出文件时使用哪一个规则!为了处理这种情况，提供了`ruleorder`解决冲突的规则

   ```
   ruleorder: rule1 > rule2 > rule3
   ```

   意思是优先使用rule1，其次是rule2，只有在前两个规则不能应用时再使用rule3

