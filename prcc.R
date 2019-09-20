#
# Copyright (C) 2019 Center for Disease Dynamics, Economics and Policy
#
# This file is part of The Hospital CRE Intervention Assessment Model (hCREiAM)
# 
# hCREiAM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# hCREiAM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with hCREiAM  If not, see <https://www.gnu.org/licenses/>.
#
#
#

# This function computes the partial rank correlation coefficient for a latin hypercube sample
# input: a dataframe

require("rms")

PRCC = function(x, sort.results = TRUE, sort.abs = TRUE) {
    x=as.data.frame(x)
    N = length(x[, 1])
    k = length(x[1, ])

    # Rank each variable\n
    r = matrix(NA, nrow = N, ncol = k)
    for (i in 1:k) {
        r[, i] = rank(x[, i])
    }

    # save output ranks
    outputvar = r[, k]

    # r is the matrix where each column (1 for each input parameter)
    # contains ranks of the simulation values for that parameter
    r = r[, -k]
    k = k - 1

    # If two of the input parameters have exactly the same ranking for
    # every run, then only one of the parameters should be used in the
    # calculation of PRCC
    dropcols = 0
    for (i in 1:(k - 1)) {
        dups = seq(1, k)
        for (j in 1:i) dups[j] = 0
        for (j in (i + 1):k) {
            a = which(r[, i] == r[, j])
            if (length(a) != N) {
                dups[j] = 0
            }
        }

        # keep track of duplicate columns
        for (j in 1:k) {
            if (dups[j] > 0) {
                dropcols = c(dropcols, j)
            }
        }
    }

    # remove 0 as a dropcol
    dropcols = dropcols[-1]

    # get unique column numbers
    dropcols = unique(dropcols)

    # sort dropcol list
    dropcols = sort(dropcols)

    # reverse list of columns to be dropped so that you
    # can drop them using the # as an index and it wont
    # interfere with subsequent drops
    dropcols = rev(dropcols)

    # drop duplicate columns
    if (length(dropcols) > 0) {
        for (i in 1:length(dropcols)) {
            r = r[, -dropcols[i]]
        }
    }

    # adjust K if you droped any columns
    k = length(r[1, ])

    # bind output rankings to input rankings matrix
    r = cbind(r, outputvar)

    ### Generate C Matrix ###
    mu = (1 + N)/2
    C.ij = matrix(NA, nrow = k + 1, ncol = k + 1)
    for (i in 1:(k + 1)) {
        for (j in 1:(k + 1)) {
            C.ij[i, j] = sum((r[, i] - mu) * (r[, j] - mu))/sqrt(sum((r[, i] -
                                                                          mu)^2) * sum((r[, j] - mu)^2))
        }
    }
    #browser()
    B = matinv(C.ij)

    gamma.ij = rep(0, k)
    t.iy = rep(0, k)
    p.iy = rep(0, k)

    for (i in 1:(k)) {
        # the PRCC between the ith input parameter and the outcome variable
        gamma.ij[i] = -B[i, k + 1]/sqrt(B[i, i] * B[k + 1, k + 1])

        # the significance of a nonzero PRCC is tested by computing t.iy
        # the distribution of t.iy approximates a students T with N-2 degrees of freedom
        t.iy[i] = gamma.ij[i] * sqrt((N - 2)/(1 - gamma.ij[i]))

        p.iy[i] = 2 * pt(-abs(t.iy[i]), df = N - 2)
    }

    # account for dropped columns
    dropcols = sort(dropcols)
    if (length(dropcols) > 0) {
        for (i in 1:length(dropcols)) {
            if (dropcols[i] > length(gamma.ij)) {
                # append to end
                gamma.ij = c(gamma.ij[1:(dropcols[i] - 1)], 0)
                t.iy = c(t.iy[1:(dropcols[i] - 1)], "dropped")
                p.iy = c(p.iy[1:(dropcols[i] - 1)], "-")
            } else {
                # insert in location
                gamma.ij = c(gamma.ij[1:(dropcols[i] - 1)], 0, gamma.ij[dropcols[i]:(length(gamma.ij))])
                t.iy = c(t.iy[1:(dropcols[i] - 1)], "dropped", t.iy[(dropcols[i]):(length(t.iy))])
                p.iy = c(p.iy[1:(dropcols[i] - 1)], "-", p.iy[(dropcols[i]):(length(p.iy))])
            }
        }
    }

    # calculate absolute value column to be used for sorting
    gamma.ij.abs = abs(gamma.ij)

    # create output dataframe
    #vals = data.frame(cbind(gamma.ij, t.iy, p.iy, gamma.ij.abs), row.names = names(x[1:(length(x[1,]) - 1)]))
    vals = data.frame(cbind(gamma.ij, t.iy, p.iy, gamma.ij.abs), row.names = colnames(x)[-length(colnames(x))])

    if (sort.results == TRUE) {
        if (sort.abs == TRUE) {
            vals = vals[order(vals$gamma.ij.abs, decreasing = TRUE), ]
        } else {
            vals = vals[order(vals$gamma.ij, decreasing = TRUE), ]
        }
    }

    # drop absolute value column
    vals = vals[, -4]

    vals

}

