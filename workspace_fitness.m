function fitness_value = workspace_fitness(arg)
%workspace_fitness 工作空间的适应度函数，假定导轨行程无限，以最小半径为适应度值
    param_errors = zeros(54,1);
    pkm = PKM(param_errors,arg(1),arg(2),arg(3),arg(4),arg(5));
    r=pkm.getWorkspaceRadius();
    fitness_value = min(r);
end

