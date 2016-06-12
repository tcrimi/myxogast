![Slime mold growth](https://tedideas.files.wordpress.com/2014/06/gif4b.gif?w=1000&h=563)

# Myxogast
## A graph aligner in Rust

Myxogast will be a network aligner based on the [classic dynamic programming](http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) [algorithms](https://en.wikipedia.org/wiki/Viterbi_algorithm), and using a novel JSON reference sequence format supporting:
 * wildcards ('N' and '*')
 * branch points
 * HMM-style probabilities

For example,
```JSON
    [{"id":    "domain_1",
     "seq":    "*ATGCATGC"},
     {"id":    "domain_2",
      "branch": [{"id":   "A",
                  "seq":  "GGCGGC"},
                 {"id":   "B",
                  "seq":  "TTATAG"}]
     },
     {"id":    "domain_3",
      "dist":  [{"A": 0.5, "T": 0.5}, {"A": 0.1, "T": 0.1, "G": 0.8}]},
     "TATTATA"]
```
