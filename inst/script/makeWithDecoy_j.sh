echo $1
export R_LIBS_SITE="/scratch/LIBWITHPROZOR/r-site-library/"
my_temp_fastaDir="tempFastaForR35"
Rscript --vanilla $R_LIBS_SITE/prozor/script/fgcz_create_fasta.R $1 --summary -o /scratch/LIBWITHPROZOR/fasta_db/tempFastaForR35/
echo "---------"
cd $my_temp_fastaDir
echo "-----looking for fasta file----"
file_name=$(ls $1_*)
echo "$file_name"

# SCP
scp "$file_name" "$USER@130.60.193.35:/srv/www/htdocs/FASTA/"
echo "doing scp to r35"
# cleaning up locally
rm "$file_name"
cd ..


#echo "$(<$1_d.txt)"
fastastore=$(ls /srv/www/htdocs/FASTA/$1_d*)
cat "$1_d.txt" | bfabric_save_fasta.py $2 $fastastore

