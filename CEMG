#!/bin/bash
####################################################################################
DATE=`date '+%Y-%m-%d-%H:%M'`
exec > >(tee -i CEMG_PIPER_${DATE}.log)
exec 2>&1
####################################################################################

Don=SHCM15
proj=$PWD
org=$PWD/raw-reads
trim=$PWD/raw-reads_trimmed
bin=$PWD/bin
assem=$PWD/assembled
binning=$PWD/assem_binned
Res=$PWD/results
human=$PWD/raw-reads_trim_Humann3
B_marker=$PWD/DonorB_markers
mkdir -p $trim $assem $Res
mkdir -p $binning $B_marker


echo "#####################################"
echo "####### Req softwares         #######"
Trimomatic="/dataone/common/software/Trimmomatic-0.38/trimmomatic-0.38.jar"
echo $Trimmomatic
Samtools="/dataone/common/software/samtools-1.11/samtools"
echo $Samtools
InterLeave=$PWD/bin/interleave_fastq.py
echo $InterLeave
Spades=$PWD/bin/SPAdes-3.15.2-Linux/bin/spades.py
echo $Spades
LongExractor=$PWD/bin/LongContigExtractor.py
echo $LongExractor
path1="/dataone/common/software/metabat2-V2021"
bin_info_sorter=$PWD/bin/bin_info_sorter.R
echo $bin_info_sorter
echo "#####################################"


echo "#######################################"
echo "Trimming the fastq files; trimmomatic"
echo "#######################################"
Adapters=$bin/NEB_Costum.fa
mkdir -p $trim
cd $org
for R1 in *_R1.fastq.gz
do
        sample=${R1%_R1*}
	R2=${sample}_R2.fastq.gz
	P1=$trim/${sample}_1P.fastq.gz
	P2=$trim/${sample}_2P.fastq.gz
	U1=$trim/${sample}_1U.fastq.gz
	U2=$trim/${sample}_2U.fastq.gz
	if [ -f $trim/${sample}_1P.fastq.gz ]; then
        echo "The file '$trim/${sample}_1P.fastq.gz' exists."
        else
        echo "The file '$trim/${sample}_1P.fastq.gz' is not found."
        java -jar $Trimomatic PE -threads 30 \
        $R1 $R2 $P1 $U1 $P2 $U2 \
        ILLUMINACLIP:$Adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
        fi
done

echo "####################################################"
echo "#######building CEMG library, co-assem pp + stool ##"
echo "####################################################"

Don_lib=$assem/${Don}_CEMG
Don_temp=$assem/${Don}_temp
mkdir -p $Don_temp

if [ -f ${Don_lib}_inter.fastq.gz ]; then
	echo "The file '${Don_lib}_inter.fastq.gz' exists."
else
        echo "The file '${Don_lib}_inter.fastq.gz' is not found."

	cd $trim
	cp *P.fastq.gz $Don_temp/
	gunzip $Don_temp/*.fastq.gz
	
	for sample_R1 in *_1P.fastq.gz; do
                sample_R2=${sample_R1%_1P*}_2P.fastq.gz
		sample=${sample_R1%_1P*}
                echo "${Don} contains: $sample_R1, $sample_R2"
                python2.7 $InterLeave $Don_temp/${sample}_1P.fastq $Don_temp/${sample}_2P.fastq $Don_temp/${sample}_inter.fastq
	done
	echo "removing SHCM15ana10b_inter.fastq from CEMG pool, this is not part of CEMG"
	mv $Don_temp/SHCM15ana10b_inter.fastq $Don_temp/SHCM15ana10b_inter_OUTofCEMG.fastq
        echo "'$Don_temp/*_inter.fastq' into $Don_lib"
        cat $Don_temp/*_inter.fastq >> ${Don_lib}_inter.fastq
        echo "Cleaning up the tem directory"
        gzip $Don_temp/*_inter.fastq
	gzip ${Don_lib}_inter.fastq
        rm $Don_temp/*_1P*.fastq
	rm $Don_temp/*_2P*.fastq

fi

echo "#####################################"
echo "#### Assembling the donor libs ######"
echo "#####################################"

if [ -f ${Don_lib}_spade/contigs.fasta ]; then
	echo "The file '${Don_lib}_spade/contigs.fasta' exists."
else
	echo "The file '${Don_lib}_spade/contigs.fasta' is not found."
	echo "assembling '${Don_lib}_inter.fastq'"
        python3 $Spades --meta --pe1-12 ${Don_lib}_inter.fastq.gz -o ${Don_lib}_spade \
                -t 150 --memory 2000
fi

echo "########################################"
echo "####### single stool assembly ## #### ##"
echo "########################################"

Don_stool=$assem/${Don}_DMG
out=${Don_stool}_spade/contigs.fasta
if [ -f $out ]; then
	echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
	python3 $Spades --meta --pe1-12 $Don_temp/${Don}Stool_inter.fastq.gz -o ${Don_stool}_spade \
        -t 50 --memory 1500
fi


echo "########################################"
echo "## Adding lib name to contig headers  ##"
echo "########################################"
assemQ=$assem/quality
mkdir -p $assemQ

out=$assemQ/${Don}_CEMG_contigs_edit_headers.txt
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."

        cd $assem
        for dir in *_spade; do
                echo $dir
                sample=${dir%_spade}
                echo $sample
                sed "/^>/s/^>/>${sample}__/g" $dir/contigs.fasta > $dir/contigs_edit.fasta
                grep -e "^>" $dir/contigs_edit.fasta > $assemQ/${sample}_contigs_edit_headers.txt
        done
fi

echo "########################################"
echo "## Comparing contigs cumulative plots ##"
echo "########################################"

out=$Res/assembly_cumulative.png
if [ -f $out ]; then
        echo "The file '$out' exists."
else
        echo "The file '$out' is not found."
        cd $assem
        Rscript $bin/CumuContigs.R quality $Res/assembly_cumulative.png
fi


echo "################################"
       	echo "Donors: $Don"
        Don_temp=$assem/${Don}_temp
        contig_file=$assem/${Don}_CEMG_spade/contigs_edit.fasta
        echo $contig_file

        echo "################################"
        echo "#####selecting large contigs####"
        echo "################################"

        dRef=${contig_file%.fasta}_1k.fa
        if [ -f $dRef ]; then
        echo "The file '$dRef' exists."
        else
        echo "The file '$dRef' is not found."
        echo "Selecting contig >= 1kb"
        python2.7 $LongExractor $contig_file
        mv ${contig_file}_out $dRef
        fi

        echo "################################"
        echo "#indexing contigs Large contigs#"
        echo "################################"

        if [ -f ${dRef}.fai ]; then
        echo "The file '${dRef}.fai' exists."
        else
        echo "The file '${dRef}.fai' is not found."
        echo "Indexing the contig file now"
        bwa index $dRef
        samtools faidx $dRef
        fi

        echo "################################"
        echo "#####Building anvio contigdb####"
        echo "################################"

	echo "changing conda env to anvio-6.2"
	#activate conda env:
	source /home/shekas3/anaconda3/bin/activate anvio-6.2
        
	dAnvi=${dRef%.fa}.db
        if [ -f $dAnvi ]; then
        echo "The file '$dAnvi' exists."
        else
        echo "The file '$dAnvi' is not found."
        echo "Building anvio contig database"
        anvi-gen-contigs-database -f $dRef -o $dAnvi
        echo "Building HMMs models"
        anvi-run-hmms -c $dAnvi
        fi

        echo "################################"
        echo "####### Short-read mapping #####"
        echo "################################"

        cover=${Don_lib}_coverage
        mkdir -p $cover
        inter=$Don_temp
	cd $inter
	for file in *_inter.fastq.gz; do
                echo "Mapping $file to $Don as follow:"
                sample_file=$inter/${file}
		sample=${file%_inter.fastq.gz}
                echo "$sample_file"
		echo $sample
	        if [ -f $cover/${sample}.bam ]; then
                echo "The file '$cover/${sample}.bam' exists."
                else
                echo "The file '$cover/${sample}.bam' is not found."
                echo "Mapping short-reads to assembly"
                bwa mem -t 15 -p $dRef $sample_file | $Samtools sort -@ 15 | $Samtools view -@ 15 -bS -o $cover/${sample}.bam
                fi

                if [ -f $cover/${sample}.bam.bai ]; then
                echo "The file '$cover/${sample}.bam.bai' exists."
                else
                echo "The file '$cover/${sample}.bam.bai' is not found."
                $Samtools index $cover/${sample}.bam
                fi

                if [ -f $cover/${sample}.coverage.txt ]; then
                echo "The file '$cover/${sample}.coverage.txt' exists."
                else
                echo "The file '$cover/${sample}.coverage.txt' is not found."
                echo "Calculating the number of mapped reads"
                $Samtools coverage $cover/${sample}.bam > $cover/${sample}.coverage.txt
                fi

                if [ -f $cover/${sample}.bam-ANVIO_PROFILE/PROFILE.db ]; then
                echo "The file '$cover/${sample}.bam-ANVIO_PROFILE/PROFILE.db' exists."
                else
                echo "The file '$cover/${sample}.bam-ANVIO_PROFILE' is not found."
                echo "Building anvio sample database"
                anvi-profile -i $cover/${sample}.bam -c $dAnvi
                fi

	done

	echo "########################################"
        echo "Binning the cotnigs using Metabat"
        echo "########################################"
	
        BIN=${binning}/${Don}_CEMG_bins
        mkdir -p $BIN

        metabat_tbl=$binning/${Don}_CEMG_depthForMetabat.txt
        if [ -f $metabat_tbl ]; then
        echo "The file '$metabat_tbl' exists."
        else
        echo "The file '$metabat_tbl' is not found."
        echo "Making a coverage table from bam file for metabat2"
        $path1/jgi_summarize_bam_contig_depths --outputDepth $metabat_tbl $cover/*.bam
        fi

        if [ -f $BIN/bin.1.fa ]; then
        echo "The file '$BIN/bin.1.fa' exists."
        else
        echo "The file '$BIN/bin.1.fa' is not found."
        echo "Binning the contigs"
        $path1/metabat2 -i $dRef -a $metabat_tbl -o $BIN/bin
        fi

        binInfo=${BIN}/contig_in_bins_anvi.txt
        if [ -f $binInfo ]; then
        echo "The file '$binInfo' exists."
        else
        echo "The file '$binInfo' is not found."
        grep -e "^>" ${BIN}/bin* > ${BIN}/contig_in_bins.txt
        Rscript $bin_info_sorter ${BIN}/contig_in_bins.txt $binInfo
        rm ${BIN}/contig_in_bins.txt
        fi

        echo "################################"
        echo "#####Building anvio sampledb####"
        echo "################################"

        dAnviProfile=$cover/SAMPLES-MERGED/PROFILE.db
        if [ -f $dAnviProfile ]; then
        echo "The file '$dAnviProfile' exists."
        else
        echo "The file '$dAnviProfile' is not found."
        anvi-merge $cover/*ANVIO_PROFILE/PROFILE.db -o $cover/SAMPLES-MERGED -c $dAnvi
        echo "Importing metabat2 bin info"
        anvi-import-collection $binInfo -p $dAnviProfile -c $dAnvi --contigs-mode -C metabat2
        anvi-script-add-default-collection -p $dAnviProfile -c $dAnvi -C default
        fi

        echo "################################"
        echo "#####Exporting coverage table###"
        echo "################################"

        dAnviSUM_M=$cover/SAMPLES-MERGED-SUM-metabat2/bins_summary.txt
        if [ -f $dAnviSUM_M ]; then
        echo "The file '$dAnviSUM_M' exists."
        else
        echo "The file '$dAnviSUM_M' is not found."
        anvi-summarize -p $dAnviProfile -c $dAnvi -C metabat2 -o $cover/SAMPLES-MERGED-SUM-metabat2 --init-gene-coverages
        fi

        dAnviSUM_D=$cover/SAMPLES-MERGED-SUM-default/bins_summary.txt
        if [ -f $dAnviSUM_D ]; then
        echo "The file '$dAnviSUM_D' exists."
        else
        echo "The file '$dAnviSUM_D' is not found."
        anvi-summarize -p $dAnviProfile -c $dAnvi -C default -o $cover/SAMPLES-MERGED-SUM-default --init-gene-coverages
        fi

        echo "#########################################"
        echo "check the quality of BINS using checkM"
        echo "##########################################"

        BINQ=${binning}/${Don}_CEMG_bins_quality
        mkdir -p $BINQ
        echo "CHANGING THE CONDA ENVIRONMENT!!!!"
        source /home/shekas3/anaconda3/bin/activate py2

        if [ -f $BINQ/lineage_wf/lineage.ms ]; then
        echo "The file '$BINQ/lineage_wf/lineage.ms' exists."
        else
        echo "The file '$BINQ/lineage_wf/lineage.ms' is not found."
        checkm lineage_wf -x fa -t 10 $BIN $BINQ/lineage_wf
        fi

        if [ -f $BINQ/lineage_wf/${Don}_bin_qc.txt ]; then
        echo "The file '$BINQ/lineage_wf/${Don}_bin_qc.txt' exists."
        else
        echo "The file '$BINQ/lineage_wf/${Don}_bin_qc.txt' is not found."
        checkm qa $BINQ/lineage_wf/lineage.ms \
        $BINQ/lineage_wf -o 1 --tab_table > $BINQ/lineage_wf/${Don}_bin_qc.txt
        fi

        echo "######################################"
        echo "Running gtdbtk (taxonomy) for each bin"
        echo "######################################"

        BINT=${binning}/${Don}_CEMG_bins_taxonomy
        mkdir -p $BINT
        echo "CHANGING THE CONDA ENVIRONMENT!!!!"
        source /home/shekas3/anaconda3/bin/activate gtdbtk

        gtdb_out="$BINT/classify/${Don}.bac120.summary.tsv"
        if [ -f $gtdb_out ]; then
        echo "The file '$gtdb_out' exists."
        else
        echo "The file '$gtdb_out' is not found."
        gtdbtk classify_wf --genome_dir $BIN --cpus 20 --out_dir $BINT --extension fa --prefix $Don
        fi


echo "#########################################################"
echo "comparing SHCM15 VS DonorB"
echo "1. percentage of mapped reads"
echo "#########################################################"

B_DMGs="/dataone/shekas3/UCFMT1/DonorB_metagenomics/Trim_interleaved"
B_cov=$assem/DonorB_DMG_coverage
mkdir -p $B_cov

out=$B_cov/DonorB_D_2013.flagstat
if [ -f $out ]; then
	echo "The file '$out' exists."
else
	echo "The file '$out' is not found."
	cd $B_DMGs
	ext="_001_inter.fastq"
	for file in *$ext
	do
		echo $file
		sample=${file%$ext}
		echo $sample
                bwa mem -t 15 -p $dRef \
		$file | $Samtools sort -@ 15 | $Samtools view -@ 15 -bS \
		-o $B_cov/${sample}.bam

		$Samtools flagstat $B_cov/${sample}.bam > $B_cov/${sample}.flagstat

	done
	cd $proj
fi


echo "#########################################################"
echo "comparing SHCM15 VS DonorB"
echo "2. Donor B's commonly engrafted strains marker coverage"
echo "#########################################################"

CEGs=$B_marker/gene_engraft_all_com_complete.fasta
SPM=$B_marker/species_markers.fa

if [ -f ${CEGs}.fai ]; then
   echo "The file '${CEGs}.fai' exists."
else
   echo "The file '${CEGs}.fai' is not found."
   echo "Indexing the contig file now"
   bwa index $CEGs
   samtools faidx $CEGs
fi

if [ -f ${SPM}.fai ]; then
   echo "The file '${SPM}.fai' exists."
else
   echo "The file '${SPM}.fai' is not found."
   echo "Indexing the contig file now"
   bwa index $SPM
   samtools faidx $SPM
fi

B_marker_cov=$B_marker/coverage
mkdir -p $B_marker_cov
stool_s=$Don_temp/${Don}Stool_inter.fastq.gz

out=$B_marker_cov/${Don}Stool_to_CEGs.coverage.txt
if [ -f $out ]; then
   echo "The file '$out' exists."
else
   echo "The file '$out' is not found."

	bwa mem -t 15 -p $CEGs $stool_s \
 	-O 60 -E 10 -L 100 | $Samtools sort -@ 15 | $Samtools view \
	-@ 15 -F 4 -o $B_marker_cov/${Don}Stool_to_CEGs.bam
	$Samtools coverage $B_marker_cov/${Don}Stool_to_CEGs.bam \
	> $B_marker_cov/${Don}Stool_to_CEGs.coverage.txt

	bwa mem -t 15 -p $SPM $stool_s \
        -O 60 -E 10 -L 100 | $Samtools sort -@ 15 | $Samtools view \
        -@ 15 -F 4 -o $B_marker_cov/${Don}Stool_to_SpeciesMarker.bam
	$Samtools coverage $B_marker_cov/${Don}Stool_to_SpeciesMarker.bam \
	> $B_marker_cov/${Don}Stool_to_SpeciesMarker.coverage.txt
fi


echo "#########################################################"
echo "Figures: comparing B vs. SHCM15"
echo "#########################################################"

out=$Res/marker_in_SHCM15.png
if [ -f $out ]; then
	echo "the file '$out' exists."
else
        echo "the file '$out' is not found."
	Rscript $bin/BvsSHCM15.R
fi





:<<"END"
echo "#########################################################"
echo "Running Humann3 on the metagenomic reads"
echo "#########################################################"
source /home/shekas3/anaconda3/bin/activate humann3
#Humann2 wants to have paired-end reads in a concatanated form, it can't properly
# use the paired-end information
mkdir -p $human

cd $inter
for sample in *_inter.fastq.gz
do
        echo $sample
        basename=${sample%_inter.fastq.gz}
        human_out=$human/${basename}/${basename}_inter_genefamilies.tsv
        if [ -f $human_out ]; then
        echo "the file '$human_out' exists."
        else
        echo "the file '$human_out' is not found."
        humann --input ${sample} --output $human/${basename} --threads 15
        fi
done


END


