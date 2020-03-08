#!/bin/bash


NO_TASKS=5
NO_TESTS=10
EXEC=tema2
VALG_FLAG=0

MAX_BONUS=10
#TEST_POINTS=(2 2 2 2 2 2 2 2 2)
TASK1_POINTS=(1 1 1 1 1 1 1 1 1 1)
TASK2_POINTS=(1 1 1 1 1 1 1 1 1 1)
TASK3_POINTS=(1 1 1 1 1 2 2 2 2 2)
TASK4_POINTS=(2 3 3 3 3 3 3 5 5 5)
TASK5_POINTS=(1 1 1 1 1 1 1 1 1 1)

#time is in ms
TASK1_TIME=(1000 1000 1000 1000 1000 1000 1000 4000 4000 4000)
TASK2_TIME=(1000 1000 1000 1000 1000 1000 1000 4000 4000 4000)
TASK3_TIME=(1000 1000 1000 1000 1000 1000 1000 4000 4000 4000)
TASK4_TIME=(3000 3000 3000 3000 3000 3000 3000 50 100 180)
TASK5_TIME=(3000 3000 3000 3000 3000 3000 3000 10000 10000 10000)
#TIMEOUT_TIME=(1 1 1 1 1 1 12 1 3 5)

VALGRIND_TIMEOUT_TIME=(14 14 14 14 14 14 14 50 50 50 150)

ERRORS_NUMBER=0
CODING_STYLE_BONUS=0
CURRENT_DIRECTORY=`pwd`
ONE=1
TWO=2
MESSAGE_SIZE=12
PASSED_TESTS_GRADE=0



echo "Compiling implementation source files..."
make || exit 1
cd _test
echo "Compiling test source files..."
make || exit 1
echo

score=""
MEMCHECK="valgrind --leak-check=full \
    --show-reachable=yes \
    --vex-iropt-register-updates=allregs-at-mem-access \
    --show-leak-kinds=all \
    --error-exitcode=100"

for i in `seq 1 7`; do
    echo "Running test $i"
    test_score="`LD_LIBRARY_PATH=.. ./test$i | tail -n1`"
    [[ $? -ne 0 ]] && test_score="0"

    echo "Running valgrind on test $i"
    LD_LIBRARY_PATH=.. $MEMCHECK "./test$i" &> /dev/null
    if [[ $? -eq 0 ]]; then
        echo "Valgrind OK"
        echo "Final score for test $i is $test_score"
        score="$score+$test_score"
    else
        echo "Valgrind FAILED"
        echo "Final score for test $i is 0"
    fi
    echo
done

score="0$score"
score="`echo "$score" | bc`"

echo "==============================================="
echo -n "Final score: "
echo "$score/80"
make clean

cd ..
make clean

function checkBonus
{   
    printf '%*s\n' "${COLUMNS:-$(($(tput cols) - $ONE))}" '' | tr ' ' -
    echo "" > checkstyle.txt
    
    echo -ne "Coding style Errors:\n"
    for entry in *.{cpp,hpp,h};
    do  
        if [[ -f $entry ]]; then
            ./cpplint.py "$entry" > checkstyle.txt
            YOUR_BONUS=`cat checkstyle.txt | grep "Total errors found: "| tr -dc '0-9'`
            ERRORS_NUMBER=$(( $ERRORS_NUMBER + $YOUR_BONUS ))
        fi
    done
    
    if [[ $ERRORS_NUMBER != "0" ]]; then
        printf '\nCoding style check: %*sFAIL' "${COLUMNS:-$(($(tput cols) - $MESSAGE_SIZE))}" '' | tr ' ' -
        echo -ne "\n\t\tYou have a total of $ERRORS_NUMBER coding style errors.\n"
    
    else
        echo -ne "\n\tNone?\n\t\tWow, next level!"
        printf '\n\nCoding style check: %*sOK' "${COLUMNS:-$(($(tput cols) - $MESSAGE_SIZE + $TWO ))}" '' | tr ' ' -
        
        CODING_STYLE_BONUS=$(echo "scale=1; $score / 8" | bc -l)
        
        # CODING_STYLE_BONUS=$(($CODING_STYLE_BONUS))
    fi 
    rm -f checkstyle.txt
}

function printBonus
{
    if [[ $ERRORS_NUMBER -gt 0 ]]; then
    printf '\n%*s' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
        echo -ne  "\n\tAnd you almost did it!\n"
    else
        echo -ne "\n\tDamn!\n"
    fi
    echo -ne "\n\t\tYou get $CODING_STYLE_BONUS/$MAX_BONUS for coding style.\n\n"
}


checkBonus
printBonus

