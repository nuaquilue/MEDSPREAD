############ Used in post.fire.r ############
count.spp <- function(x){
  return(c(sum(x==1), sum(x==2), sum(x==3), sum(x==4), sum(x==5), sum(x==6), sum(x==7),
           sum(x==8), sum(x==9)))
}

select.cohort <- function(x){
  if(sum(x)==0)
    return(10)
  else
    return(sample(1:10, 1, replace=F, prob=x))
}

