import numpy as np
import Poisson_template as program
# Run the code implemented in (b) to solve the two test problems defined in (a). 
# Solve for the mesh resolutions N∈[20,40,80,160,320,640,1280] 
# and compute the maximum norm of the error 
# maxi |ui−u(xi)|. What convergence order do you get?


# 1️. 테스트 문제 선택
testproblems = ["const", "sin"]  # 두 가지 문제
parameters = program.define_default_parameters()

# 2️. N values given
N_values = [20, 40, 80, 160, 320, 640, 1280]

for problem_name in testproblems:
    print(f"\nSolving test problem: {problem_name}")
    testproblem = program.get_testproblem(problem_name)
    
    err_max_list = []

    for N in N_values:
        # solver 실행
        # err_max, _ = program.my_driver("full", testproblem, parameters, N) # full
        err_max, _ = program.my_driver("sparse", testproblem, parameters, N)

        err_max_list.append(err_max)
    
    # 최대오차 출력
    print("N\tMax error")
    for N, err in zip(N_values, err_max_list):
        print(f"{N}\t{err:.2e}")
    
    # 수렴율 계산 (log-log)
    err_max_array = np.array(err_max_list)
    rate = np.log(err_max_array[:-1]/err_max_array[1:])/np.log(2)
    print("Convergence rates:", rate)
