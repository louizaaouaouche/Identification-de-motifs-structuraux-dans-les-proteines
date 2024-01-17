filename=$1
aff=$2

var1=".fasta"
var2=".pdb"
var3="rcsb_pdb_"

echo "Dimension 1:"
up=$(echo $filename | tr '[a-z]' '[A-Z]')
fd1=$var3$up$var1
echo $fd1
python3 main.py $fd1 1 $2


echo -e "_____________________________\n"


echo "Dimension 2:"
low=$(echo $filename | tr '[A-Z]' '[a-z]')
fd2=$low$var2
echo $fd2
python3 main.py $fd2 2 $2
