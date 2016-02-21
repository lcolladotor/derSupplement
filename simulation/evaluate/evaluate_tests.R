foo <- function(ex) {
    ex <- unlist(trans)
    trans_id <- rep(names(trans), elementLengths(trans))
    idx <- match(trans_id, chosen$ucsckg_id)
    ex$DEr1 <- chosen$rep1[idx] != 'normal'
    ex$DEr2 <- chosen$rep2[idx] != 'normal'
    ex$DEr3 <- chosen$rep3[idx] != 'normal'
    
    ex$case1 <- chosen$case1[idx]
    ex$case2 <- chosen$case2[idx]
    ex$case3 <- chosen$case3[idx]
    
    return(ex)
}

exons_foo <- foo()

percDE(exons_foo)

tables_exons_f <- mapply(count_comp, stats_GR, replicates, types,
    SIMPLIFY = FALSE, MoreArgs = list(reference = exons_foo, cut = 0.05))
empirical_exons_f <- emp(tables_exons_f)
index_exons_f <- mapply(index_comp, stats_GR, replicates, types,
    SIMPLIFY = FALSE, MoreArgs = list(reference = exons_foo, cut = 0.05))
case_result_exons_f <- mapply(case_result, index_exons_f, replicates,
    SIMPLIFY = FALSE, MoreArgs = list(reference = exons_foo))
emp_sum(empirical_exons_f)


exons_u <- processExons(unique(unlist(trans[elementLengths(trans) < 3])))
exons_u <- exons_disjoin[countOverlaps(exons_disjoin, trans) == 1]
exons_u <- exons[countOverlaps(exons, trans) == 1]
percDE(exons_u)

tables_exons_u <- mapply(count_comp, stats_GR, replicates, types,
    SIMPLIFY = FALSE, MoreArgs = list(reference = exons_u, cut = 0.05))
empirical_exons_u <- emp(tables_exons_u)
index_exons_u <- mapply(index_comp, stats_GR, replicates, types,
    SIMPLIFY = FALSE, MoreArgs = list(reference = exons_u, cut = 0.05))
case_result_exons_u <- mapply(case_result, index_exons_u, replicates,
    SIMPLIFY = FALSE, MoreArgs = list(reference = exons_u))
emp_sum(empirical_exons_u)

case_result_exons_u[['regionMatrix.1.edgeR']]
