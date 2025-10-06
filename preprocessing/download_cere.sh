#!/usr/bin/env bash
#SBATCH --job-name=sample_python_job    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=2gb                     # Job memory request (per node)
#SBATCH --time=78:05:00               # Time limit hrs:min:sec
#SBATCH --output=rclone_test_%j.log   # Standard output and error log
#SBATCH --partition celltypes         # Partition used for processing
#SBATCH --tmp=10G                     # Request the amount of space your jobs needs on /scratch/fast

cd ~/Matthew/data/taxtest/extra/zhong_cerebellum_PRJNA695270/fastqs/

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/067/SRR13565267/SRR13565267_1.fastq.gz -o SRR13565267_GSM5047776_GW16-01_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/067/SRR13565267/SRR13565267_2.fastq.gz -o SRR13565267_GSM5047776_GW16-01_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/072/SRR13565272/SRR13565272_1.fastq.gz -o SRR13565272_GSM5047781_GW21-02_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/072/SRR13565272/SRR13565272_2.fastq.gz -o SRR13565272_GSM5047781_GW21-02_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/069/SRR13565269/SRR13565269_1.fastq.gz -o SRR13565269_GSM5047778_GW18-01_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/069/SRR13565269/SRR13565269_2.fastq.gz -o SRR13565269_GSM5047778_GW18-01_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/073/SRR13565273/SRR13565273_1.fastq.gz -o SRR13565273_GSM5047782_GW21-03_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/073/SRR13565273/SRR13565273_2.fastq.gz -o SRR13565273_GSM5047782_GW21-03_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/071/SRR13565271/SRR13565271_1.fastq.gz -o SRR13565271_GSM5047780_GW21-01_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/071/SRR13565271/SRR13565271_2.fastq.gz -o SRR13565271_GSM5047780_GW21-01_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/070/SRR13565270/SRR13565270_1.fastq.gz -o SRR13565270_GSM5047779_GW18-02_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/070/SRR13565270/SRR13565270_2.fastq.gz -o SRR13565270_GSM5047779_GW18-02_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/066/SRR13565266/SRR13565266_1.fastq.gz -o SRR13565266_GSM5047775_GW14-02_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/066/SRR13565266/SRR13565266_2.fastq.gz -o SRR13565266_GSM5047775_GW14-02_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/064/SRR13565264/SRR13565264_1.fastq.gz -o SRR13565264_GSM5047773_GW12-02_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/064/SRR13565264/SRR13565264_2.fastq.gz -o SRR13565264_GSM5047773_GW12-02_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/063/SRR13565263/SRR13565263_1.fastq.gz -o SRR13565263_GSM5047772_GW12-01_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/063/SRR13565263/SRR13565263_2.fastq.gz -o SRR13565263_GSM5047772_GW12-01_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/068/SRR13565268/SRR13565268_1.fastq.gz -o SRR13565268_GSM5047777_GW16-02_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/068/SRR13565268/SRR13565268_2.fastq.gz -o SRR13565268_GSM5047777_GW16-02_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/065/SRR13565265/SRR13565265_1.fastq.gz -o SRR13565265_GSM5047774_GW14-01_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/065/SRR13565265/SRR13565265_2.fastq.gz -o SRR13565265_GSM5047774_GW14-01_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/076/SRR13565276/SRR13565276_1.fastq.gz -o SRR13565276_GSM5047785_GW27_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/076/SRR13565276/SRR13565276_2.fastq.gz -o SRR13565276_GSM5047785_GW27_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/074/SRR13565274/SRR13565274_1.fastq.gz -o SRR13565274_GSM5047783_GW21-04_cerebellum_Homo_sapiens_RNA-Seq_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR135/074/SRR13565274/SRR13565274_2.fastq.gz -o SRR13565274_GSM5047783_GW21-04_cerebellum_Homo_sapiens_RNA-Seq_2.fastq.gz
