print_green () {
  tput setaf 2
  echo $@
  tput setaf 7
}
export -f print_green

print_red () {
  tput setaf 1
  echo $@
  tput setaf 7
}
export -f print_red

echo "Our smaller match interval runs from SNP 6 to SNP 189, so we test -bits in increments of 20 up to 81."
echo "Germline finds the homoz match tract at -bits up to (including) 94. At bits 95, the first slice is SNPs 1-95, the second slice is SNPs 96-190, and the third slice is SNPs 191-200."
echo "Since none of these slices are fully homozygous, Germline doesn't find the homoz tract to extend from."
export BITS_TO_TEST=(11 21 41 61 81)
