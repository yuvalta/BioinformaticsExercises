import re


def viterbi_seq():
    viterbi_output = []
    prob_dict = {'red': 0, 'yellow': 0, 'green': 0}

    # with open("viterbi_seq.txt") as fileObj:
    #
    #     for word in fileObj:
    #         viterbi_output.append(word)

    f = open('viterbi_seq.txt', 'r')
    for word in f:
        viterbi_output = re.split(',', word)

    for state in range(len(viterbi_output)):
        prob_dict[viterbi_output[state]] = prob_dict[viterbi_output[state]] + 1

    print(prob_dict)
    # {'red': 2967, 'yellow': 1999, 'green': 5035}


def calc_probs_of_random_letters():
    prob_dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    with open("random_seq.txt") as fileObj:
        for line in fileObj:
            for ch in line:
                prob_dict[ch] = prob_dict[ch] + 1

    print(prob_dict)
    # {'A': 2537, 'C': 2542, 'G': 2432, 'T': 2489}

# calc_probs_of_random_letters()
# viterbi_seq()
