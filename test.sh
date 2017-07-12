if [ "$#" -ne 1 ]; then
    echo -e "\e[101m Input number of processes for testing \e[0m"
    exit
fi

for script in lagrange_primal energy lagrange_mixed
do
    echo -e "\e[44m Testing "$script" \e[0m";
    mpirun -np $1 python $script.py;
done

# Graph dofs vs iters
python plot_convergence.py
