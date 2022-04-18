set -euxo pipefail

all_tests=($(find ../test/out/graph500*))

for test in ${all_tests[@]}; do
    echo -e "\n\033[0;33m RUN: ${test} \033[0m\n"
    ./${test} || exit 1
done