function result = splicemats(ID, data)
% 输入参数：
% ID: 1*N cell类型，表示组学数据的ID
% data: 20530*N double类型，表示组学数据
% filename: string类型，表示输出csv文件的名称
% 输出参数：无
% 转置ID和组学数据
ID = cellfun(@char, ID, 'UniformOutput', false); %这一行将元胞元胞数组转换成字符串元胞数组
ID = ID'; % 转置元胞数组，使行向量转变成列向量
ID = table(ID) % 将元胞转换成table，这是关键步骤！
data = data'; % 将data转置
data = table(data) % 将data转换成table，这更是关键步骤！
% 拼接ID和组学数据，只有两个table可以不同数据类型拼接在一起！
result = [ID, data];
end



% 记住，使用writetable(D, 'met_ID.csv')这种格式来写入文件！
% 生存数据单独的处理方法：
patient_ID = table(patient_ID)
patient_dmfs_e = table(patient_dmfs_e)
patient_dmfs_time = table(patient_dmfs_time)
Sur = [patient_ID,patient_dmfs_e,patient_dmfs_time]