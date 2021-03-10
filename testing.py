read = "CTACAGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGTCTCAGCGTCAGTCACTGTCCAGAAGGCCGCCTTCGCCACCGGTGTTCCGCTTGATATCTACGAATTTCACCTCTACACCAAGCATTCCACCTTCCTCTCCAGGACTCTAGCCCTGCAGTTTCGAGTGCACTTCCACGGTTGAGCCGTGGGATTTCACACCCGACTTGCAGGGCCGCCTACACGCGCTTTACGCCCAGTAATTCCGAACAACGTTTGCACCCTCTGTCTTACCGCGGCTGCTGGCAC"
seq = "GCCGGTTGCTTGATCCACCAGGAACCGTCTGTTAGCACCGCCGAACGCCTCGTTTTCCCCCTGCTTCTGCACGTTTGAGACACGCGTGAACGCCCTCGCCCAAGCCGCTACCCCCTTCCCAAGTGCTGGAATTCTCTGGGCTCGGATGCTGGATGCCATCCGCCAGGAGAAGTTCGACTACGCGCTGAGATGGCTGGAGCGGATGCGCCCGCTGGAGGTGCGTGACGGCGCCCTGGTGATGGGCGTTCCGGACCGCTTCTTCCGCGACTGGGTGGATGACCACTACCGCCCCATGCTGGATGTCCAGCTCGCGCGGATGGGGGAGGGCCTGACCTCCATCGCCTATGAGGTGGTGGAGGGCCTGGAGCCGGACCCGCATTTTCCACCCACACCCTCAGTCAAGGCGAGCGCCACGCGCCCGGGCCGGCTGAACGCGCGCTTCACCTTCGACACCTTCGTGGTGGCGGACAGCAACCAGCTCCCGGCCGCCGCGGCGCAGGCCGTGGCCGACAAGCCGGGCCACCACTACAACCCGCTCTACATCTACGGCGGCACGGGCCTGGGCAAGACGCACCTGCTCCAGGCGGTGGGCAACCTCATCTGGGAGCGGGATCCGTCCCAGCGCGTGGTGTTCCTCTCCAGCGAGCAGTTCACAAACGAGTACGTGGAGAGCGTCCGCGAGCACCGCATGGGGGACTTCCGCCGGAAGTTCCGTGAGGAGTGCGACGTGCTCCTCCTCGACGACATCCAGTTCCTCGGCAAGCGTGAGGAGACGCAGAAGGAGTTCTTCTACACCTTCAACACGCTCTACGAGATGAACAAGGCCATCGTCCTCACCAGCGACACCGTGCCCGCGGAGATTCCGGGCCTGGAGGACCGGCTGCGCAGCCGCTTCGCCATGGGGCTGATGACGGACATCCGCGAGCCCACCTACGAGACGCGGGTGGCCATCCTCCAGAAAAAGGCCGTGGCCGAAGGTCTGGACCTCCCGGACTCGGTGGCGCACTTCATCGCCAAGCACATCCAGAAGAACGTACGCGAGCTGGAAGGCGCGCTGGTGAAGCTGTCCGCGGTGCACAGCCTGACGCGCCAGCCGGTGACGGAGGATTTCGCGTCCCAGGTGCTGCGCGACATCCTCCCCGCTCACAGCGCGGTGGACGTGGAGTCCATCCAGCGCGAGGTGGCCCGCTTCTACAAAGTCACGGTGGAGTCGCTGAAGGAAGACCGCCGGCACAAGGCCCTGGCCCATGCGCGACAGGTGGCCATGTACCTCAGCCGCAAGCTGACGAAGAGCTCCTTCCCGGAAATCGCCGCGCGCTTCAGCAAGGACCACTCCACCGTCATCTCCGCGGTGCGCAAGGTGGAAGGCCTCCGCATGACGGACGCCACCGTGCAGCGCGAGTTGGCGGAGCTGGAATTGAAGCTCGGCAATCCCTGAAGCAGCCGTGCCTTTGTTTCATGTCGGACGA"


from sw_algorithm import sw_cpp

print(sw_cpp(read="AGGTTCGACGTT", seq="AGGTCGGTCAA", algorithm="sw", seq_type="nuc"))

print(sw_cpp(read="GCCCAGTGGT", seq="AGTTTGCCCAGTGGTGCCCACA", algorithm="sw", seq_type="aa")) 

print(sw_cpp(read="GCCCAGTGGT", seq="AGTTTGCCCAGTGGTGCCCACA", algorithm="nw", seq_type="aa"))

print(sw_cpp(read="GCCCAGTGGT", seq="AGTTTGCCCAGTGGTGCCCACA", algorithm="nw", seq_type="nuc"))


