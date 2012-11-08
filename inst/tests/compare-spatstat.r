context('spatstat equivalence for r, K, L')

test_that('bronzefilter gives same r and K for both implementations', {
    pattern <- spatstat::bronzefilter

    our_res <- summ_func(pattern, edge_correction = 'translational',
                         mtf_name = '1', do_besags_L = TRUE)
    our_a <- our_res[['a']]
    our_K <- our_a['original', 'K_1', , drop=TRUE]
    our_L <- our_a['original', 'L_1', , drop=TRUE]
    our_r <- our_res[['r']]

    sp_res <- spatstat::Kest(pattern)
    sp_K <- sp_res[['trans']]
    sp_L <- spatstat::Lest(pattern)[['trans']]
    sp_r <- sp_res[['r']]

    expect_that(our_K, equals(sp_K))
    expect_that(our_L, equals(sp_L))
    expect_that(our_r, equals(sp_r))
})
