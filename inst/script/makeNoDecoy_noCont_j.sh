echo $1
export R_LIBS_SITE="/scratch/LIBWITHPROZOR/r-site-library/"
my_temp_fastaDir="tempFastaForR35"
Rscript --vanilla $R_LIBS_SITE/prozor/script/fgcz_create_fasta.R $1 nodecoy --summary -c none -o /scratch/LIBWITHPROZOR/fasta_db/tempFastaForR35/
echo "-----looking for fasta file----"
echo "$my_temp_fastaDir"
#file_name=$(find "$my_temp_fastaDir" -type f -name "*.fasta")
cd $my_temp_fastaDir
file_name=$(ls *.fasta)
echo "$file_name"
#echo "I am $USER"

scp "$file_name" "$USER@130.60.193.35:/srv/www/htdocs/FASTA/"
echo "doing scp to r35"
# cleaning up locally
rm "$file_name"
cd ..


#echo "$(<$1_d.txt)"
#fastastore=$(ls /srv/www/htdocs/FASTA/$1*)
fastastore=$(find /srv/www/htdocs/FASTA/ -type f -name "$1_*")
cat "$1.txt" | bfabric_save_fasta.py $2 $fastastore

