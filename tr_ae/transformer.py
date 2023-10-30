import torch.nn as nn
from datasets import *

# the latent dimension
d_feature = 3

# the dimension of model or the size of the data segment.
d_model = 128

# the number of heads of Multi-Head Attention
n_heads = 8

# the dimension of inner-layer of Position-wise Feed-Forward Networks
d_ff = 2048

"""
The input  "Scaled Dot-Product Attention"
consists of queries and keys of dimension d_k, and values of dimension d_v.
"""
d_k = d_v = int(d_model / n_heads)

# The encoder/decoder is composed of a stack of N identical layers.
n_layers = 6

# Dimensions of each cell's genetic data
cell_dim = 512


# Set the parameters of the model
def set_tr_parameter(d_model_num, train_set, latent_dim):
    global d_model
    d_model = d_model_num
    global d_k
    d_k = int(d_model / n_heads)
    global d_v
    d_v = int(d_model / n_heads)
    global cell_dim
    cell_dim = train_set.shape[1] * train_set.shape[2]
    global d_feature
    d_feature = latent_dim


# Positional Encoding
class PositionalEncoding(nn.Module):
    def __init__(self, d_model, dropout=0.1, max_len=5000):
        super(PositionalEncoding, self).__init__()
        self.dropout = nn.Dropout(p=dropout)
        pos_table = np.array([
            [pos / np.power(10000, 2 * i / d_model) for i in range(d_model)]
            if pos != 0 else np.zeros(d_model) for pos in range(max_len)])
        pos_table[1:, 0::2] = np.sin(pos_table[1:, 0::2])
        pos_table[1:, 1::2] = np.cos(pos_table[1:, 1::2])
        self.pos_table = torch.FloatTensor(pos_table).cuda()

    def forward(self, enc_inputs):
        enc_inputs += self.pos_table[:enc_inputs.size(1), :]
        return self.dropout(enc_inputs.cuda())


# Scaled Dot-Product Attention
class ScaledDotProductAttention(nn.Module):
    def __init__(self):
        super(ScaledDotProductAttention, self).__init__()

    def forward(self, Q, K, V):
        scores = torch.matmul(Q, K.transpose(-1, -2)) / np.sqrt(d_k)
        attn = nn.Softmax(dim=-1)(scores)
        context = torch.matmul(attn, V)
        return context


# Multi-Head Attention
class MultiHeadAttention(nn.Module):
    def __init__(self):
        super(MultiHeadAttention, self).__init__()
        self.W_Q = nn.Linear(d_model, d_k * n_heads, bias=False)
        self.W_K = nn.Linear(d_model, d_k * n_heads, bias=False)
        self.W_V = nn.Linear(d_model, d_v * n_heads, bias=False)
        self.fc = nn.Linear(n_heads * d_v, d_model, bias=False)

    def forward(self, input_Q, input_K, input_V):
        residual, batch_size = input_Q, input_Q.size(0)
        # n_heads = 8
        Q = self.W_Q(input_Q).view(batch_size, -1, n_heads, d_k).transpose(1, 2)  # Q: []
        K = self.W_K(input_K).view(batch_size, -1, n_heads, d_k).transpose(1, 2)  # K: []
        V = self.W_V(input_V).view(batch_size, -1, n_heads, d_v).transpose(1, 2)  # V: []

        context = ScaledDotProductAttention()(Q, K, V)
        context = context.transpose(1, 2).reshape(batch_size, -1, n_heads * d_v)
        output = self.fc(context)
        return nn.LayerNorm(d_model).cuda()(output + residual)


#  Position-wise Feed-Forward Networks
class PoswiseFeedForwardNet(nn.Module):
    def __init__(self):
        super(PoswiseFeedForwardNet, self).__init__()
        self.fc = nn.Sequential(
            nn.Linear(d_model, d_ff, bias=False),
            nn.ReLU(),
            nn.Linear(d_ff, d_model, bias=False))

    def forward(self, inputs):
        residual = inputs
        output = self.fc(inputs)
        return nn.LayerNorm(d_model).cuda()(output + residual)


class EncoderLayer(nn.Module):
    def __init__(self):
        super(EncoderLayer, self).__init__()
        self.enc_self_attn = MultiHeadAttention()
        self.pos_ffn = PoswiseFeedForwardNet()

    def forward(self, enc_inputs):
        enc_outputs = self.enc_self_attn(enc_inputs, enc_inputs, enc_inputs)
        enc_outputs = self.pos_ffn(enc_outputs)
        return enc_outputs


class Encoder(nn.Module):
    def __init__(self):
        super(Encoder, self).__init__()
        self.en_embedding = nn.Linear(d_model, d_model, bias=False)
        self.pos_emb = PositionalEncoding(d_model)
        self.layers = nn.ModuleList([EncoderLayer() for _ in range(n_layers)])

        # Reduce dimension:from cell_dim to latent dimension
        self.dec_d = nn.Linear(cell_dim, d_feature, bias=False)

    def forward(self, enc_inputs):
        enc_outputs = self.pos_emb(enc_inputs)
        for layer in self.layers:
            enc_outputs = layer(enc_outputs)

        batch_size = enc_outputs.size(0)
        enc_outputs = enc_outputs.reshape(batch_size, -1)
        enc_outputs = self.dec_d(enc_outputs)
        return enc_outputs


class DecoderLayer(nn.Module):
    def __init__(self):
        super(DecoderLayer, self).__init__()
        self.dec_enc_attn = MultiHeadAttention()
        self.pos_ffn = PoswiseFeedForwardNet()

    def forward(self, dec_inputs):
        dec_outputs = self.dec_enc_attn(dec_inputs, dec_inputs, dec_inputs)
        dec_outputs = self.pos_ffn(dec_outputs)
        return dec_outputs


class Decoder(nn.Module):
    def __init__(self):
        super(Decoder, self).__init__()

        # Ascending dimension:from latent dimension to cell_dim
        self.fc_de = nn.Sequential(
            nn.Linear(d_feature, cell_dim, bias=False),
            nn.ReLU()
        )
        self.layers = nn.ModuleList([DecoderLayer() for _ in range(n_layers)])

    def forward(self, enc_outputs):
        enc_outputs = self.fc_de(enc_outputs)

        batch_size = enc_outputs.size(0)
        dec_outputs = enc_outputs.reshape(batch_size, -1, d_model)

        for layer in self.layers:
            dec_outputs = layer(dec_outputs)

        dec_outputs = torch.exp(dec_outputs)
        return dec_outputs


class Transformer(nn.Module):
    """
    This class implements a Transformer autoencoder
    """
    def __init__(self):
        super(Transformer, self).__init__()
        self.Encoder = Encoder().cuda()
        self.Decoder = Decoder().cuda()

    def forward(self, enc_inputs):
        enc_outputs = self.Encoder(enc_inputs)
        dec_outputs = self.Decoder(enc_outputs)
        return dec_outputs
