# Andy's nice-ish y axis labels with ggplot
# Converts 1e6 into 10^6 etc., and avoids things like 10^0 or .5 x 10^-1
# Defines a function called within scale_y_continuous, as the argument to 'labels'.
# Example:
# library(ggplot2)
# p=ggplot() + xlim(0,1)+ scale_y_continuous(limits=c(8e-4, 2e5),  trans="log10", breaks = as.vector(c(1,5) %o% 10^(-8:8)),
# minor_breaks = as.vector((1:9) %o% 10^(-8:8)), labels = ggplot_scientific_notation_axes_labels)

ggplot_scientific_notation_axes_labels <- function(l) {
  # R function 'prettyNum' puts things into scientific notation when digit strings get past a certain length,
  # but keeps things  >0.0001 and < 10^7 in non-sci notation
  l <- prettyNum(l)
  # … but for axis labels, we want to switch to powers of 10 for anything >=10^3 or < 0.01
  # (though you decide - edit on line 17)
  for(i in 1:length(l)){
    number=eval(parse(text=l)[i])
    if(!is.na(number) & (abs(number)<0.01 | abs(number)>=1000)) l[i]=format(number, scientific=TRUE)
  }
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (e^+2 -> e^2)
  l <- gsub("e\\+","e",l)  
  l <- gsub("e0","e",l)  
  l <- gsub("e-0","e-",l) 
  # turn the 'e' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  # convert 0.000.. * 10^x to 0
  l <- gsub("\\'0[\\.0]*\\'\\%\\*\\%", "", l)
  l <- gsub("10\\^0", "0", l)
  # return this as an expression
  parse(text=l)
}


library(ggplot2)
p=ggplot() + xlim(0,1) + theme_bw() + 
  scale_y_continuous(limits=c(8e-4, 2e5), trans="log10",breaks = as.vector(c(1,5) %o% 10^(-8:8)),
                     minor_breaks = as.vector((1:9) %o% 10^(-8:8)), labels = ggplot_scientific_notation_axes_labels)

#print(p)


