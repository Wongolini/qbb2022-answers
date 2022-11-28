array1=( "vietnam" "Germany" "Argentina" )
array2=( "Asia" "Europe" "America")

for i in ${!array1[@]}; do
	echo ${array1[i]} ${array2[i]}
done
