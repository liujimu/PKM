function fitness_value = workspace_fitness(arg)
%workspace_fitness �����ռ����Ӧ�Ⱥ������ٶ������г����ޣ�����С�뾶Ϊ��Ӧ��ֵ
    param_errors = zeros(54,1);
    pkm = PKM(param_errors,arg(1),arg(2),arg(3),arg(4),arg(5));
    r=pkm.getWorkspaceRadius();
    fitness_value = min(r);
end

