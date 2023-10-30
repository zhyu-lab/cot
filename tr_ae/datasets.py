import torch
import torch.utils.data as Data
import numpy as np


# The data is segmented according to "d_seg", and the part that is not divided is 0
def get_seg_data(data_ori, d_seg):
    data_ori = torch.from_numpy(data_ori)
    tensorData_split = data_ori
    train_set = None
    for i in range(data_ori.shape[0]):
        row, tensorData_split = tensorData_split.split([1, tensorData_split.shape[0] - 1], dim=0)
        total = None
        for j in range(int(data_ori.shape[1] // d_seg)):
            a1, row = row.split([d_seg, row.shape[1] - d_seg], dim=1)
            if total is None:
                total = a1
            else:
                total = torch.cat([total, a1], dim=0)

        if row.shape[1] != 0:
            # Fill the data with zeros
            temp = torch.zeros(1, d_seg - row.shape[1])
            row = torch.cat([row, temp], dim=1)
            total = torch.cat([total, row], dim=0)

        total = torch.stack([total], dim=0)
        if train_set is None:
            train_set = total
        else:
            train_set = torch.cat([train_set, total], dim=0)
    return train_set


# Cell gene data set
class CellDataSet(Data.Dataset):
    def __init__(self, enc_inputs):
        super(CellDataSet, self).__init__()
        self.enc_inputs = enc_inputs

    def __len__(self):
        return self.enc_inputs.shape[0]

    def __getitem__(self, idx):
        return self.enc_inputs[idx]
