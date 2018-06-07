# DO some shit for trajectory choices

def connect(thing_list, mode):
    this = ''
    if mode == 'start':
        for thing in thing_list:
            this += (thing)

    if mode == 'end':
        for thing in thing_list:
            this += thing
            # if thing == thing_list[-1]:
            #     this += thing
            # else:
            #     this += ('-' + thing)
    return this


# def get_solutions(A,B):
#
#     solutions = []
#     for start in A:
#         for end in B:
#             start_connect = connect(start, 'start')
#             end_connect = connect(end, 'end')
#             solutions.append(start_connect+end_connect)
#     return solutions

def get_solutions(A, B):
    solutions = []
    for start in A:
        for end in B:
            start_connect = connect(start, 'start')
            end_connect = connect(end, 'end')
            solutions.append(start_connect + end_connect)
    return solutions


def print_things(solution):
    for i in solution:
        print(i)


E = 'E'
V = 'V'
M = 'M'
J = 'J'
S = 'S'
P = 'P'

_1Ai = [[E, V],
        [E, V, E],
        [E, V, V],
        [E, V, E, V],
        [E, V, E, E],
        [E, V, V, E],
        [E, V, V, V],
        [E, M],
        [E, M, E],
        [E, M, M],
        [E, M, E, M],
        [E, M, E, E],
        [E, M, M, E],
        [E, M, M, M],
        [E, V, M],
        [E, M, V],
        [E, V, M, V],
        [E, V, M, M],
        [E, V, V, M],
        [E, M, V, M],
        [E, M, V, V],
        [E, M, M, V]]

_1B = [[J, P],
       [S, P],
       [J, S, P],
       [P]]

solutions_1 = get_solutions(_1Ai, _1B)

#######################################################################################################

_2A = [[E, E],
       [E, V],
       [E, M]]

_2B = [[J, P],
       [S, P],
       [J, S, P]]

solutions_2 = get_solutions(_2A, _2B)
#######################################################################################################

_3 = [[E, J, P],
      [E, J, S, P],
      [E, S, P]]
_3B = ['']
solutions_3 = get_solutions(_3, _3B)

# solutions = get_solutions(_1Aiaa, _1B)
# for solution in solutions:
#     print(solution)

# print_things(solutions_1)
# print_things(solutions_2)
# print_things(solutions_3)

all_solutions_dups = solutions_1 + solutions_2 + solutions_3
all_solutions = set(all_solutions_dups)
all_sols_refined_1 = []
for sol in all_solutions:
    if ((V + V + V) not in sol) and ((M + M + M) not in sol) and ((E + E + E) not in sol):
        if not (sol[0] == E and sol[1] == E):
            if (not ((J or S) in sol) and (len(sol) < 5)) or (((J or S) in sol) and (len(sol) < 6)):
                if ((V + V) not in sol) and ((M + M) not in sol) and ((E + E) not in sol):
                    if ((J + S) not in sol):
                        all_sols_refined_1.append(sol)

# print(len(all_sols_refined_1))

sols2send = []
for sol in all_sols_refined_1:
    newsol = []
    for bod in sol:
        if bod == E: newsol.append('earth')
        if bod == V: newsol.append('venus')
        if bod == M: newsol.append('mars')
        if bod == J: newsol.append('jupiter')
        if bod == S: newsol.append('saturn')
        if bod == P: newsol.append('pluto')
    sols2send.append(newsol)


def split_list(a_list):
    half = len(a_list) / 2
    half = int(half)
    # print(half)
    return a_list[:half], a_list[half:]


sols2send.sort(key=lambda x: ''.join(x))
sols2send_half1, sols2send_half2 = split_list(sols2send)

sols2send_1, sols2send_2 = split_list(sols2send_half1)
sols2send_3, sols2send_4 = split_list(sols2send_half2)

#
print(sols2send_1)
print(sols2send_2)
print(sols2send_3)
print(sols2send_4)
