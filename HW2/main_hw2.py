import random


def create_random_data():
    dna_letters = ["A", "C", "G", "T"]
    seq_list = []

    f = open("./random_seq.txt", "w")

    for i in range(10000):
        letter = random.choice(dna_letters)
        f.write(letter)
        seq_list.append(letter)

    f.close()
    return seq_list

def viterbi(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    for i in states:
        V[0][i] = start_p[i] * emit_p[i][obs[0]]
    # Run Viterbi when t > 0
    for t in range(1, len(obs)):
        V.append({})
        for y in states:
            (prob, state) = max((V[t - 1][y0] * trans_p[y0][y] * emit_p[y][obs[t]], y0) for y0 in states)
            V[t][y] = prob
        # for i in dptable(V):
        #     print(i)
        opt = []
        for j in V:
            for x, y in j.items():
                if j[x] == max(j.values()):
                    opt.append(x)
    # the highest probability
    h = max(V[-1].values())
    print('The steps of states are -- ' + ' '.join(opt) + ' -- with highest probability of %s' % h)

    f = open("./viterbi_seq.txt", "w")
    f.write('The steps of states are -- ' + ' '.join(opt))
    f.close()


def dptable(V):
    yield " ".join(("%10d" % i) for i in range(len(V)))
    for y in V[0]:
        yield "%.7s: " % y+" ".join("%.7s" % ("%f" % v[y]) for v in V)

# random_list = create_random_data() # generate data only once

random_list = []
viterbi_output = []

with open("./random_seq.txt") as fileObj:
    for line in fileObj:
        for ch in line:
            random_list.append(ch)

states = ["green", "yellow", "red"]
start_p = {"green": 1 / 3, "yellow": 1 / 3, "red": 1 / 3}
observations = ("A", "C", "G", "T")
trans_p = {
    "green": {"green": 0.3, "yellow": 0.6, "red": 0.1},
    "yellow": {"green": 0.7, "yellow": 0.1, "red": 0.2},
    "red": {"green": 0.3, "yellow": 0.2, "red": 0.5}
}
emit_p = {
    "green": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
    "yellow": {"A": 0.7, "C": 0.1, "G": 0.1, "T": 0.1},
    "red": {"A": 0.1, "C": 0.2, "G": 0.3, "T": 0.4}
}

viterbi(tuple(random_list), states, start_p, trans_p, emit_p)



