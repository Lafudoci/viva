# Virus Integrity and Variant Analyzer Report
## Meta
Task name : test

Task ID : test_run_202104200748

Task start time : 2021-04-20 15:48:16 UTC+8

Task finish time : 2021-04-20 15:50:47 UTC+8
## Input Reads
| Reads | File name |
| ----- | --------- |
| R1 | /home/leftc/viva/test_data/AdV_R1.fastq.gz |
| R2 | /home/leftc/viva/test_data/AdV_R2.fastq.gz |


| Reads | MD5 hash |
| ----- | --------- |
| R1 | 31a20d72323904847d40e480d68143d2 |
| R2 | 59383681f63b5a478054340275632c9c |
## Input Reference
Reference by user : Yes

Origin file path : /home/leftc/viva/test_data/adv_multi_ref.fasta

Numeber of reference : 3

SPAdes assembly mode : N/A

Best contig name : N/A

Best contig identity : N/A %
| Order | Header | Length |
| ----- | ------ | ------ |
| 1 | AC_000008.1 Human adenovirus 5, complete genome | 35938 |
| 2 | M73260.1 Mastadenovirus h5 gene, complete genome | 35935 |
| 3 | KF268127.1 Human adenovirus C strain human/USA/CL_42/1988/5[P5H5F5], complete genome | 35931 |
## Reads Filter Statistics
### Filter
|   | Before | After |
| - | ------ | ----- |
| Total reads | 155844 | 155266 |
| Total base  | 22870134 | 22787703 |
| Q30 rate    | 0.943333 | 0.944953 |
| Mean length R1 | 146 | 146 |
| Mean length R2 | 146 | 146 |
### Duplication
Duplication rate : 2.34%
## Host Genome Removal
Remove genome: N/A

Percentage of removed reads: N/A
## Alignment
### Mapping rate

#### Reference #1 :AC_000008.1 Human adenovirus 5, complete genome
| Aligner | Overall mapped rate |
| ------- | ------------------- |
| Bowtie2 | 99.56% |
| BWA MEM | 99.55% |
#### Reference #2 :M73260.1 Mastadenovirus h5 gene, complete genome
| Aligner | Overall mapped rate |
| ------- | ------------------- |
| Bowtie2 | 99.56% |
| BWA MEM | 99.55% |
#### Reference #3 :KF268127.1 Human adenovirus C strain human/USA/CL_42/1988/5[P5H5F5], complete genome
| Aligner | Overall mapped rate |
| ------- | ------------------- |
| Bowtie2 | 99.56% |
| BWA MEM | 99.55% |
### Coverage

#### Reference #1 :AC_000008.1 Human adenovirus 5, complete genome
| Aligner | Start base | End base | Covered base | Mean depth |
| ------- | ---------- | -------- | ------------ | ---------- |
| Bowtie2 | 1 | 35938 | 99.9666 % | 630.471 x |
| BWA MEM | 1 | 35938 | 99.9666 % | 630.175 x |
#### Reference #2 :M73260.1 Mastadenovirus h5 gene, complete genome
| Aligner | Start base | End base | Covered base | Mean depth |
| ------- | ---------- | -------- | ------------ | ---------- |
| Bowtie2 | 1 | 35935 | 99.9666 % | 630.42 x |
| BWA MEM | 1 | 35935 | 99.9666 % | 630.133 x |
#### Reference #3 :KF268127.1 Human adenovirus C strain human/USA/CL_42/1988/5[P5H5F5], complete genome
| Aligner | Start base | End base | Covered base | Mean depth |
| ------- | ---------- | -------- | ------------ | ---------- |
| Bowtie2 | 1 | 35931 | 99.9666 % | 630.667 x |
| BWA MEM | 1 | 35931 | 99.9666 % | 630.333 x |
## Variant Calling
### Reference #1 :AC_000008.1 Human adenovirus 5, complete genome
#### LoFreq
| Position | REF | ALT | Bowtie2(Filter/AF/DP/GQ) | BWA(Filter/AF/DP/GQ) |
| -------- | --- | --- | ------- | --- |
| 4952 | G | C | - / - / - / - | PASS / 96.57% / 554 / 19897 |
| 8783 | G | A | - / - / - / - | PASS / 93.59% / 577 / 19391 |
| 11284 | T | C | - / - / - / - | PASS / 97.13% / 522 / 18217 |
| 14073 | T | TA | - / - / - / - | PASS / 20.16% / 620 / 1784 |
| 17387 | G | C | - / - / - / - | PASS / 97.44% / 507 / 18009 |
| 18754 | GA | G | - / - / - / - | PASS / 92.42% / 594 / 20765 |
| 19483 | T | A | - / - / - / - | PASS / 98.04% / 764 / 28763 |
| 19513 | T | A | - / - / - / - | PASS / 98.95% / 765 / 28636 |
| 19657 | G | A | - / - / - / - | PASS / 97.76% / 760 / 27730 |
| 19658 | A | G | - / - / - / - | PASS / 98.18% / 768 / 28017 |
| 20378 | T | C | - / - / - / - | PASS / 97.09% / 755 / 27442 |
| 21163 | C | T | - / - / - / - | PASS / 97.71% / 656 / 24284 |
| 21630 | G | A | - / - / - / - | PASS / 96.21% / 580 / 20982 |
| 25995 | A | T | - / - / - / - | PASS / 88.11% / 597 / 19169 |
| 26726 | GGCGGCA | G | - / - / - / - | PASS / 90.76% / 476 / 17262 |
| 27161 | C | T | - / - / - / - | PASS / 96.96% / 624 / 23022 |
| 27314 | C | A | - / - / - / - | PASS / 97.32% / 597 / 21351 |
| 27339 | T | C | - / - / - / - | PASS / 97.98% / 593 / 22238 |
| 27650 | T | C | - / - / - / - | PASS / 95.67% / 669 / 24058 |
| 27651 | C | T | - / - / - / - | PASS / 95.37% / 669 / 23769 |
| 28120 | T | C | - / - / - / - | PASS / 97.38% / 686 / 25095 |
| 29823 | C | CT | - / - / - / - | PASS / 30.26% / 780 / 6362 |
| 30301 | A | C | - / - / - / - | PASS / 94.98% / 617 / 21641 |
| 30302 | C | G | - / - / - / - | PASS / 93.23% / 606 / 20799 |
| 30303 | G | C | - / - / - / - | PASS / 94.03% / 620 / 21532 |
| 30403 | C | A | - / - / - / - | PASS / 95.18% / 643 / 22556 |
| 30404 | A | C | - / - / - / - | PASS / 97.11% / 658 / 23842 |
| 34343 | G | GT | - / - / - / - | PASS / 17.27% / 718 / 1808 |
| 35776 | A | C | - / - / - / - | PASS / 95.35% / 752 / 26498 |

#### Varscan2
| Position | REF | ALT | Bowtie2(Filter/AF/DP/GQ) | BWA(Filter/AF/DP/GQ) |
| -------- | --- | --- | ------- | --- |
| 4952 | G | C | - / - / - / - | PASS / 100% / 477 / 255 |
| 8783 | G | A | - / - / - / - | PASS / 100% / 482 / 255 |
| 11284 | T | C | - / - / - / - | PASS / 100% / 439 / 255 |
| 14073 | T | TA | - / - / - / - | PASS / 20.87% / 516 / 255 |
| 17387 | G | C | - / - / - / - | PASS / 100% / 426 / 255 |
| 18754 | GA | G | - / - / - / - | PASS / 91% / 489 / 255 |
| 19483 | T | A | - / - / - / - | PASS / 99.85% / 662 / 255 |
| 19513 | T | A | - / - / - / - | PASS / 99.85% / 663 / 255 |
| 19657 | G | A | - / - / - / - | PASS / 100% / 664 / 255 |
| 19658 | A | G | - / - / - / - | PASS / 100% / 671 / 255 |
| 20378 | T | C | - / - / - / - | PASS / 100% / 633 / 255 |
| 21163 | C | T | - / - / - / - | PASS / 100% / 553 / 255 |
| 21630 | G | A | - / - / - / - | PASS / 100% / 490 / 255 |
| 25995 | A | T | - / - / - / - | PASS / 100% / 499 / 255 |
| 26726 | GGCGGCA | G | - / - / - / - | PASS / 90.62% / 408 / 255 |
| 27161 | C | T | - / - / - / - | PASS / 100% / 532 / 255 |
| 27314 | C | A | - / - / - / - | PASS / 100% / 492 / 255 |
| 27339 | T | C | - / - / - / - | PASS / 100% / 487 / 255 |
| 27650 | T | C | - / - / - / - | PASS / 100% / 547 / 255 |
| 27651 | C | T | - / - / - / - | PASS / 100% / 549 / 255 |
| 28120 | T | C | - / - / - / - | PASS / 100% / 587 / 255 |
| 29823 | C | CT | - / - / - / - | PASS / 29.16% / 613 / 255 |
| 30301 | A | C | - / - / - / - | PASS / 96.41% / 534 / 255 |
| 30302 | C | G | - / - / - / - | PASS / 96.28% / 523 / 255 |
| 30303 | G | C | - / - / - / - | PASS / 96.38% / 536 / 255 |
| 30403 | C | A | - / - / - / - | PASS / 99.43% / 549 / 255 |
| 30404 | A | C | - / - / - / - | PASS / 99.46% / 562 / 255 |
| 35776 | A | C | - / - / - / - | PASS / 100% / 548 / 255 |
### Reference #2 :M73260.1 Mastadenovirus h5 gene, complete genome
#### LoFreq
| Position | REF | ALT | Bowtie2(Filter/AF/DP/GQ) | BWA(Filter/AF/DP/GQ) |
| -------- | --- | --- | ------- | --- |
| 4952 | G | C | - / - / - / - | PASS / 96.57% / 554 / 19897 |
| 8783 | G | A | - / - / - / - | PASS / 93.59% / 577 / 19391 |
| 11284 | T | C | - / - / - / - | PASS / 97.13% / 522 / 18217 |
| 14073 | T | TA | - / - / - / - | PASS / 20.16% / 620 / 1784 |
| 17387 | G | C | - / - / - / - | PASS / 97.44% / 507 / 18009 |
| 18754 | GA | G | - / - / - / - | PASS / 92.42% / 594 / 20765 |
| 19483 | T | A | - / - / - / - | PASS / 98.04% / 764 / 28763 |
| 19513 | T | A | - / - / - / - | PASS / 98.95% / 765 / 28636 |
| 19657 | G | A | - / - / - / - | PASS / 97.76% / 760 / 27730 |
| 19658 | A | G | - / - / - / - | PASS / 98.18% / 768 / 28017 |
| 20378 | T | C | - / - / - / - | PASS / 97.09% / 755 / 27442 |
| 21163 | C | T | - / - / - / - | PASS / 97.71% / 656 / 24284 |
| 21630 | G | A | - / - / - / - | PASS / 96.21% / 580 / 20982 |
| 25995 | A | T | - / - / - / - | PASS / 88.11% / 597 / 19169 |
| 26726 | GGCGGCA | G | - / - / - / - | PASS / 90.76% / 476 / 17262 |
| 27161 | C | T | - / - / - / - | PASS / 96.96% / 624 / 23022 |
| 27314 | C | A | - / - / - / - | PASS / 97.32% / 597 / 21351 |
| 27339 | T | C | - / - / - / - | PASS / 97.98% / 593 / 22238 |
| 27650 | T | C | - / - / - / - | PASS / 95.67% / 669 / 24058 |
| 27651 | C | T | - / - / - / - | PASS / 95.37% / 669 / 23769 |
| 28120 | T | C | - / - / - / - | PASS / 97.38% / 686 / 25095 |
| 29823 | C | CT | - / - / - / - | PASS / 30.26% / 780 / 6362 |
| 30301 | A | C | - / - / - / - | PASS / 94.98% / 617 / 21641 |
| 30302 | C | G | - / - / - / - | PASS / 93.23% / 606 / 20799 |
| 30303 | G | C | - / - / - / - | PASS / 94.03% / 620 / 21532 |
| 30403 | C | A | - / - / - / - | PASS / 95.18% / 643 / 22556 |
| 30404 | A | C | - / - / - / - | PASS / 97.11% / 658 / 23842 |
| 34343 | G | GT | - / - / - / - | PASS / 17.27% / 718 / 1808 |
| 34860 | T | A | - / - / - / - | PASS / 97.54% / 650 / 24278 |
| 34933 | CT | C | - / - / - / - | PASS / 97.61% / 753 / 29077 |
| 35320 | T | TCA | - / - / - / - | PASS / 91.83% / 710 / 25668 |
| 35509 | T | TC | - / - / - / - | PASS / 95.28% / 848 / 29740 |
| 35522 | C | CA | - / - / - / - | PASS / 100.00% / 825 / 29413 |
| 35773 | A | C | - / - / - / - | PASS / 95.35% / 752 / 26498 |

#### Varscan2
| Position | REF | ALT | Bowtie2(Filter/AF/DP/GQ) | BWA(Filter/AF/DP/GQ) |
| -------- | --- | --- | ------- | --- |
| 4952 | G | C | - / - / - / - | PASS / 100% / 477 / 255 |
| 8783 | G | A | - / - / - / - | PASS / 100% / 482 / 255 |
| 11284 | T | C | - / - / - / - | PASS / 100% / 439 / 255 |
| 14073 | T | TA | - / - / - / - | PASS / 20.87% / 516 / 255 |
| 17387 | G | C | - / - / - / - | PASS / 100% / 426 / 255 |
| 18754 | GA | G | - / - / - / - | PASS / 91% / 489 / 255 |
| 19483 | T | A | - / - / - / - | PASS / 99.85% / 662 / 255 |
| 19513 | T | A | - / - / - / - | PASS / 99.85% / 663 / 255 |
| 19657 | G | A | - / - / - / - | PASS / 100% / 664 / 255 |
| 19658 | A | G | - / - / - / - | PASS / 100% / 671 / 255 |
| 20378 | T | C | - / - / - / - | PASS / 100% / 633 / 255 |
| 21163 | C | T | - / - / - / - | PASS / 100% / 553 / 255 |
| 21630 | G | A | - / - / - / - | PASS / 100% / 490 / 255 |
| 25995 | A | T | - / - / - / - | PASS / 100% / 499 / 255 |
| 26726 | GGCGGCA | G | - / - / - / - | PASS / 90.62% / 408 / 255 |
| 27161 | C | T | - / - / - / - | PASS / 100% / 532 / 255 |
| 27314 | C | A | - / - / - / - | PASS / 100% / 492 / 255 |
| 27339 | T | C | - / - / - / - | PASS / 100% / 487 / 255 |
| 27650 | T | C | - / - / - / - | PASS / 100% / 547 / 255 |
| 27651 | C | T | - / - / - / - | PASS / 100% / 549 / 255 |
| 28120 | T | C | - / - / - / - | PASS / 100% / 587 / 255 |
| 29823 | C | CT | - / - / - / - | PASS / 29.16% / 613 / 255 |
| 30301 | A | C | - / - / - / - | PASS / 96.41% / 534 / 255 |
| 30302 | C | G | - / - / - / - | PASS / 96.28% / 523 / 255 |
| 30303 | G | C | - / - / - / - | PASS / 96.38% / 536 / 255 |
| 30403 | C | A | - / - / - / - | PASS / 99.43% / 549 / 255 |
| 30404 | A | C | - / - / - / - | PASS / 99.46% / 562 / 255 |
| 34860 | T | A | - / - / - / - | PASS / 100% / 559 / 255 |
| 34933 | CT | C | - / - / - / - | PASS / 96.94% / 623 / 255 |
| 34934 | T | C | - / - / - / - | PASS / 53.33% / 633 / 60 |
| 35320 | T | TCA | - / - / - / - | PASS / 91.37% / 592 / 255 |
| 35509 | T | TC | - / - / - / - | PASS / 93.63% / 660 / 255 |
| 35522 | C | CA | - / - / - / - | PASS / 96.74% / 830 / 255 |
| 35773 | A | C | - / - / - / - | PASS / 100% / 548 / 255 |
### Reference #3 :KF268127.1 Human adenovirus C strain human/USA/CL_42/1988/5[P5H5F5], complete genome
#### LoFreq
| Position | REF | ALT | Bowtie2(Filter/AF/DP/GQ) | BWA(Filter/AF/DP/GQ) |
| -------- | --- | --- | ------- | --- |
| 14073 | T | TA | - / - / - / - | PASS / 20.16% / 620 / 1784 |
| 29816 | C | CT | - / - / - / - | PASS / 30.26% / 780 / 6362 |
| 34336 | G | GT | - / - / - / - | PASS / 17.27% / 718 / 1808 |

#### Varscan2
| Position | REF | ALT | Bowtie2(Filter/AF/DP/GQ) | BWA(Filter/AF/DP/GQ) |
| -------- | --- | --- | ------- | --- |
| 14073 | T | TA | - / - / - / - | PASS / 20.87% / 516 / 255 |
| 29816 | C | CT | - / - / - / - | PASS / 29.16% / 613 / 255 |

## Draft Genome

### Reference #1 :AC_000008.1 Human adenovirus 5, complete genome

FASTA was saved to : /home/leftc/viva/tasks/test_run_202104200748/draft_genome/test_run_202104200748_draft_1.fasta

Apllied SNV : G4952C, G8783A, T11284C, G17387C, GA18754G, T19483A, T19513A, G19657A, A19658G, T20378C, C21163T, G21630A, A25995T, GGCGGCA26726G, C27161T, C27314A, T27339C, T27650C, C27651T, T28120C, A30301C, C30302G, G30303C, C30403A, A30404C, A35776C

Conflict calling : None

Mismatch calling : None
### Reference #2 :M73260.1 Mastadenovirus h5 gene, complete genome

FASTA was saved to : /home/leftc/viva/tasks/test_run_202104200748/draft_genome/test_run_202104200748_draft_2.fasta

Apllied SNV : G4952C, G8783A, T11284C, G17387C, GA18754G, T19483A, T19513A, G19657A, A19658G, T20378C, C21163T, G21630A, A25995T, GGCGGCA26726G, C27161T, C27314A, T27339C, T27650C, C27651T, T28120C, A30301C, C30302G, G30303C, C30403A, A30404C, T34860A, CT34933C, T35320TCA, T35509TC, C35522CA, A35773C

Conflict calling : None

Mismatch calling : None
### Reference #3 :KF268127.1 Human adenovirus C strain human/USA/CL_42/1988/5[P5H5F5], complete genome

FASTA was saved to : /home/leftc/viva/tasks/test_run_202104200748/draft_genome/test_run_202104200748_draft_3.fasta

Apllied SNV : None

Conflict calling : None

Mismatch calling : None
## Commands
| Duration(s) | Executed command |
| ---- | ------- |
| 2 | fastp -i /home/leftc/viva/tasks/test_run_202104200748/reads/original/test_run_202104200748_R1.fastq.gz -I /home/leftc/viva/tasks/test_run_202104200748/reads/original/test_run_202104200748_R2.fastq.gz -o /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R1.fastq.gz -O /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R2.fastq.gz -j /home/leftc/viva/tasks/test_run_202104200748/reads/fastp.json -h /home/leftc/viva/tasks/test_run_202104200748/reads/fastp.html -f 0 -t 0 -F 0 -T 0 -w 6 |
| 0 | bowtie2-build --threads 6 test_run_202104200748_ref_1.fasta test_run_202104200748_ref_1 |
| 0 | bowtie2-build --threads 6 test_run_202104200748_ref_2.fasta test_run_202104200748_ref_2 |
| 0 | bowtie2-build --threads 6 test_run_202104200748_ref_3.fasta test_run_202104200748_ref_3 |
| 11 | bowtie2 -x /home/leftc/viva/tasks/test_run_202104200748/alignment/bowtie2/test_run_202104200748_ref_1 -1 /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R1.fastq.gz -2 /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R2.fastq.gz -S test_run_202104200748_ref_1.sam -p 6 --very-sensitive-local --un-conc-gz test_run_202104200748_ref_1_unmapped_R%.fastq.gz |
| 12 | bowtie2 -x /home/leftc/viva/tasks/test_run_202104200748/alignment/bowtie2/test_run_202104200748_ref_2 -1 /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R1.fastq.gz -2 /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R2.fastq.gz -S test_run_202104200748_ref_2.sam -p 6 --very-sensitive-local --un-conc-gz test_run_202104200748_ref_2_unmapped_R%.fastq.gz |
| 11 | bowtie2 -x /home/leftc/viva/tasks/test_run_202104200748/alignment/bowtie2/test_run_202104200748_ref_3 -1 /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R1.fastq.gz -2 /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R2.fastq.gz -S test_run_202104200748_ref_3.sam -p 6 --very-sensitive-local --un-conc-gz test_run_202104200748_ref_3_unmapped_R%.fastq.gz |
| 1 | samtools sort -@ 6 test_run_202104200748_ref_1.sam -o test_run_202104200748_ref_1.sorted.bam |
| 0 | samtools index -@ 6 test_run_202104200748_ref_1.sorted.bam |
| 0 | samtools sort -@ 6 test_run_202104200748_ref_2.sam -o test_run_202104200748_ref_2.sorted.bam |
| 0 | samtools index -@ 6 test_run_202104200748_ref_2.sorted.bam |
| 1 | samtools sort -@ 6 test_run_202104200748_ref_3.sam -o test_run_202104200748_ref_3.sorted.bam |
| 0 | samtools index -@ 6 test_run_202104200748_ref_3.sorted.bam |
| 0 | bwa index -p test_run_202104200748_ref_1 test_run_202104200748_ref_1.fasta |
| 0 | bwa index -p test_run_202104200748_ref_2 test_run_202104200748_ref_2.fasta |
| 0 | bwa index -p test_run_202104200748_ref_3 test_run_202104200748_ref_3.fasta |
| 1 | bwa mem -t 6 /home/leftc/viva/tasks/test_run_202104200748/alignment/bwa/test_run_202104200748_ref_1 /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R1.fastq.gz /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R2.fastq.gz -o test_run_202104200748_ref_1.sam |
| 2 | bwa mem -t 6 /home/leftc/viva/tasks/test_run_202104200748/alignment/bwa/test_run_202104200748_ref_2 /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R1.fastq.gz /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R2.fastq.gz -o test_run_202104200748_ref_2.sam |
| 2 | bwa mem -t 6 /home/leftc/viva/tasks/test_run_202104200748/alignment/bwa/test_run_202104200748_ref_3 /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R1.fastq.gz /home/leftc/viva/tasks/test_run_202104200748/reads/test_run_202104200748_R2.fastq.gz -o test_run_202104200748_ref_3.sam |
| 0 | samtools sort -@ 6 test_run_202104200748_ref_1.sam -o test_run_202104200748_ref_1.sorted.bam |
| 0 | samtools index -@ 6 test_run_202104200748_ref_1.sorted.bam |
| 1 | samtools sort -@ 6 test_run_202104200748_ref_2.sam -o test_run_202104200748_ref_2.sorted.bam |
| 0 | samtools index -@ 6 test_run_202104200748_ref_2.sorted.bam |
| 0 | samtools sort -@ 6 test_run_202104200748_ref_3.sam -o test_run_202104200748_ref_3.sorted.bam |
| 0 | samtools index -@ 6 test_run_202104200748_ref_3.sorted.bam |
| 0 | samtools flagstat -@ 6 test_run_202104200748_ref_1.sorted.bam |
| 0 | samtools flagstat -@ 6 test_run_202104200748_ref_2.sorted.bam |
| 0 | samtools flagstat -@ 6 test_run_202104200748_ref_3.sorted.bam |
| 0 | samtools flagstat -@ 6 test_run_202104200748_ref_1.sorted.bam |
| 0 | samtools flagstat -@ 6 test_run_202104200748_ref_2.sorted.bam |
| 1 | samtools flagstat -@ 6 test_run_202104200748_ref_3.sorted.bam |
| 0 | samtools coverage test_run_202104200748_ref_1.sorted.bam |
| 1 | samtools coverage test_run_202104200748_ref_2.sorted.bam |
| 0 | samtools coverage test_run_202104200748_ref_3.sorted.bam |
| 1 | samtools coverage test_run_202104200748_ref_1.sorted.bam |
| 0 | samtools coverage test_run_202104200748_ref_2.sorted.bam |
| 1 | samtools coverage test_run_202104200748_ref_3.sorted.bam |
| 2 | lofreq faidx test_run_202104200748_ref_1.fasta |
| 6 | lofreq call-parallel --pp-threads 6 -f test_run_202104200748_ref_1.fasta -o test_run_202104200748_bowtie2_ref_1_lofreq.vcf --call-indels -N -B -q 20 -Q 20 -m 20 test_run_202104200748_ref_1.indelqual.sorted.bam |
| 2 | lofreq faidx test_run_202104200748_ref_2.fasta |
| 6 | lofreq call-parallel --pp-threads 6 -f test_run_202104200748_ref_2.fasta -o test_run_202104200748_bowtie2_ref_2_lofreq.vcf --call-indels -N -B -q 20 -Q 20 -m 20 test_run_202104200748_ref_2.indelqual.sorted.bam |
| 2 | lofreq faidx test_run_202104200748_ref_3.fasta |
| 6 | lofreq call-parallel --pp-threads 6 -f test_run_202104200748_ref_3.fasta -o test_run_202104200748_bowtie2_ref_3_lofreq.vcf --call-indels -N -B -q 20 -Q 20 -m 20 test_run_202104200748_ref_3.indelqual.sorted.bam |
| 2 | lofreq faidx test_run_202104200748_ref_1.fasta |
| 6 | lofreq call-parallel --pp-threads 6 -f test_run_202104200748_ref_1.fasta -o test_run_202104200748_bwa_ref_1_lofreq.vcf --call-indels -N -B -q 20 -Q 20 -m 20 test_run_202104200748_ref_1.indelqual.sorted.bam |
| 2 | lofreq faidx test_run_202104200748_ref_2.fasta |
| 5 | lofreq call-parallel --pp-threads 6 -f test_run_202104200748_ref_2.fasta -o test_run_202104200748_bwa_ref_2_lofreq.vcf --call-indels -N -B -q 20 -Q 20 -m 20 test_run_202104200748_ref_2.indelqual.sorted.bam |
| 2 | lofreq faidx test_run_202104200748_ref_3.fasta |
| 6 | lofreq call-parallel --pp-threads 6 -f test_run_202104200748_ref_3.fasta -o test_run_202104200748_bwa_ref_3_lofreq.vcf --call-indels -N -B -q 20 -Q 20 -m 20 test_run_202104200748_ref_3.indelqual.sorted.bam |
| 1 | samtools mpileup -B -f /home/leftc/viva/tasks/test_run_202104200748/alignment/bowtie2/test_run_202104200748_ref_1.fasta /home/leftc/viva/tasks/test_run_202104200748/alignment/bowtie2/test_run_202104200748_ref_1.sorted.bam |
| 8 | java -jar /home/leftc/bioapp/varscan2/VarScan.v2.4.4.jar mpileup2cns --min-avg-qual 20 --P-value 0.01 --output-vcf 1 |
| 2 | samtools mpileup -B -f /home/leftc/viva/tasks/test_run_202104200748/alignment/bowtie2/test_run_202104200748_ref_2.fasta /home/leftc/viva/tasks/test_run_202104200748/alignment/bowtie2/test_run_202104200748_ref_2.sorted.bam |
| 9 | java -jar /home/leftc/bioapp/varscan2/VarScan.v2.4.4.jar mpileup2cns --min-avg-qual 20 --P-value 0.01 --output-vcf 1 |
| 1 | samtools mpileup -B -f /home/leftc/viva/tasks/test_run_202104200748/alignment/bowtie2/test_run_202104200748_ref_3.fasta /home/leftc/viva/tasks/test_run_202104200748/alignment/bowtie2/test_run_202104200748_ref_3.sorted.bam |
| 7 | java -jar /home/leftc/bioapp/varscan2/VarScan.v2.4.4.jar mpileup2cns --min-avg-qual 20 --P-value 0.01 --output-vcf 1 |
| 2 | samtools mpileup -B -f /home/leftc/viva/tasks/test_run_202104200748/alignment/bwa/test_run_202104200748_ref_1.fasta /home/leftc/viva/tasks/test_run_202104200748/alignment/bwa/test_run_202104200748_ref_1.sorted.bam |
| 8 | java -jar /home/leftc/bioapp/varscan2/VarScan.v2.4.4.jar mpileup2cns --min-avg-qual 20 --P-value 0.01 --output-vcf 1 |
| 2 | samtools mpileup -B -f /home/leftc/viva/tasks/test_run_202104200748/alignment/bwa/test_run_202104200748_ref_2.fasta /home/leftc/viva/tasks/test_run_202104200748/alignment/bwa/test_run_202104200748_ref_2.sorted.bam |
| 8 | java -jar /home/leftc/bioapp/varscan2/VarScan.v2.4.4.jar mpileup2cns --min-avg-qual 20 --P-value 0.01 --output-vcf 1 |
| 1 | samtools mpileup -B -f /home/leftc/viva/tasks/test_run_202104200748/alignment/bwa/test_run_202104200748_ref_3.fasta /home/leftc/viva/tasks/test_run_202104200748/alignment/bwa/test_run_202104200748_ref_3.sorted.bam |
| 7 | java -jar /home/leftc/bioapp/varscan2/VarScan.v2.4.4.jar mpileup2cns --min-avg-qual 20 --P-value 0.01 --output-vcf 1 |

## Versions
| Tool | Version |
| ---- | ------- |
| fastp | 0.21.0 |
| samtools | 1.11 |
| bowtie2 | 2.4.2 |
| bwa | 0.7.17-r1198-dirty |
| lofreq | 2.1.5 |
| varscan2 | 2.4.4 |
| spades | v3.15.0 |
| last_commit | v0.7.0-50-g7408c27 |
