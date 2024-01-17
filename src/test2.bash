filename=$1

var1=".fasta"
var2=".pdb"
var3="rcsb_pdb_"

echo "Dimension 2:"
low=$(echo $filename | tr '[A-Z]' '[a-z]')
fd2=$low$var2
echo $fd2
python3 main.py $fd2 2 $2