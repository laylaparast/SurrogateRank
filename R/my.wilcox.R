#in delong paper x = a, y = b
my.wilcox = function(x,y) {
  if(x>y) return(1)
  if(x==y) return(1/2)
  if(x<y) return(0)
}