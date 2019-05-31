for i in *_out
do
        sed -i '1,5d' $i
        sort -n -k 1,1 -o $i $i
done
