function result = splicemats(ID, data)
% ���������
% ID: 1*N cell���ͣ���ʾ��ѧ���ݵ�ID
% data: 20530*N double���ͣ���ʾ��ѧ����
% filename: string���ͣ���ʾ���csv�ļ�������
% �����������
% ת��ID����ѧ����
ID = cellfun(@char, ID, 'UniformOutput', false); %��һ�н�Ԫ��Ԫ������ת�����ַ���Ԫ������
ID = ID'; % ת��Ԫ�����飬ʹ������ת���������
ID = table(ID) % ��Ԫ��ת����table�����ǹؼ����裡
data = data'; % ��dataת��
data = table(data) % ��dataת����table������ǹؼ����裡
% ƴ��ID����ѧ���ݣ�ֻ������table���Բ�ͬ��������ƴ����һ��
result = [ID, data];
end



% ��ס��ʹ��writetable(D, 'met_ID.csv')���ָ�ʽ��д���ļ���
% �������ݵ����Ĵ�������
patient_ID = table(patient_ID)
patient_dmfs_e = table(patient_dmfs_e)
patient_dmfs_time = table(patient_dmfs_time)
Sur = [patient_ID,patient_dmfs_e,patient_dmfs_time]