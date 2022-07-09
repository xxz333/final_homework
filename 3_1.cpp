#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<algorithm>
#include<vector>
#include<arm_neon.h>
using namespace std;
//矩阵行数
const int row_num_1 = 22;
const int row_num_2 = 8;
int row_n_1 = row_num_1;//消元子
int row_n_2 = row_num_2;//被消元行
int col_num = 130;//矩阵列数
//用于存放读入的十进制的消元子和被消元行
int** a, ** b;
void malloc_ab();
//将数据读入内存
void read_data(string file_name, int** x);
//将数据转化为二进制形式存储
unsigned int** bin_a;
unsigned int** bin_b;
int bin_col_num;//存储二进制数的矩阵列数   如：将一行的130个数分为32+32+32+32+2的5个4Byte的形式来存储，则bin_col_num=5;
void malloc_ab_bin();
//将十进制数转为二进制存储 eg.5->100000，哪一位为1就代表存在与该位位数相同的十进制数的值
void dec_to_bin(int** x, unsigned int** bin_x, int row_num);
//将二进制数转为十进制
void bin_to_dec(int** x, vector<unsigned int*> bin_v);
//存储二进制转化为十进制的结果
int** a1, ** b1;
void malloc_a1b1();
//将申请的内存们释放
void d_mem();
//将消元子和被消元行保存到vector中
vector<unsigned int*> bin_a_v;
vector<unsigned int*> bin_b_v;
//找到一行中首个非零元素位置(bin_a/bin_b)
int find_first(unsigned int* k);
//进行两行间的消去操作
void handle_one(unsigned int* k1, unsigned int* k2, int n_bei, bool& judge1, bool& judge2);
//查看被消元行首项是否有与该位置对应的消元子
bool equal_loc(int first_1_loc);
//找到被消元行插入到消元子中时的插入位置
int find_insert_loc(int& loc);
bool cmp1(unsigned int* a1, unsigned int* a2)
{
    if (a1[0] > a2[0])
        return true;
    else
        return false;
}

int main()
{
    //读取文件中的数据
    string str1 = "消元子.txt";
    string str2 = "被消元行.txt";
    malloc_ab();
    read_data(str1, a);
    read_data(str2, b);
    //将数据转化为二进制形式存储
    malloc_ab_bin();
    dec_to_bin(a, bin_a, row_num_1);
    dec_to_bin(b, bin_b, row_num_2);
    //将消元子和被消元行保存到vector中
    bin_a_v.insert(bin_a_v.begin(), bin_a, bin_a + row_num_1);
    bin_b_v.insert(bin_b_v.begin(), bin_b, bin_b + row_num_2);
    //对于每个被消元行，遍历所有消元子,第一次消元结束之后继续遍历
    for (int i = 0; i < bin_b_v.size(); i++)
        //for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < bin_a_v.size(); j++)
            //for (int j = 0; j < 1; j++)
        {
            //在将原来的被消元行变为消元子时出错
            int n_1 = find_first(bin_b_v[i]);
            int n_2 = find_first(bin_a_v[j]);
            //判断此时消元子与被消元行的首元素是否位于同一位置
            if (n_1 == n_2 && n_1 != -1)//位于同一位置
            {
                bool judge1 = false, judge2 = false;
                handle_one(bin_a_v[j], bin_b_v[i], i, judge1, judge2);//两行间进行消去
                if (judge1 == true)//可以直接跳到下一行
                {
                    break;
                }
                if (judge2 == true)
                    break;
            }
        }
    }
    //将二进制数据转为十进制形式存储
    malloc_a1b1();
    bin_to_dec(b1, bin_b_v);
    bin_to_dec(a1, bin_a_v);
    cout << "被消元行" << endl;
    for (int i = 0; i < bin_b_v.size(); i++)
    {
        for (int j = 0;; j++)
        {
            if (b1[i][j] == -1)
                break;
            else
                cout << b1[i][j] << " ";
        }
        cout << endl;
    }
    //内存释放
    d_mem();
}

void malloc_ab()
{
    //分配空间存储原始数据
    a = new int* [row_num_1];
    b = new int* [row_num_2];
    for (int i = 0; i < row_num_1; i++)
        a[i] = new int[col_num];
    for (int i = 0; i < row_num_2; i++)
        b[i] = new int[col_num];
}
void read_data(string file_name, int** x)
{
    ifstream ifs(file_name, ifstream::in);
    string str;
    int row_num_now = 0;
    while (getline(ifs, str))
    {
        int intS[3000];
        //int i = 0;//某个数的位数
        int num = 0;//一行中的数据数目
        int len = str.length();
        int len_temp = 0;
        string temp[3000];
        for (; len_temp < len; )
        {
            while (str[len_temp] != ' ' && len_temp < len)
            {
                temp[num] += str[len_temp];
                len_temp++;
            }
            //将string类型的temp转化为int类型
            intS[num] = atoi(temp[num].c_str());
            num++;
            len_temp++;
        }
        for (int i = 0; i < num; i++)
            x[row_num_now][i] = intS[i];
        x[row_num_now][num] = -1;//结束符
        row_num_now++;
    }
    ifs.close();
}
void malloc_ab_bin()
{
    //分配空间
    bin_a = new unsigned int* [row_num_1];
    bin_b = new unsigned int* [row_num_2];
    if ((col_num % 32) == 0)
        bin_col_num = col_num / 32;
    else
        bin_col_num = col_num / 32 + 1;
    for (int i = 0; i < row_num_1; i++)
        bin_a[i] = new unsigned int[bin_col_num];
    for (int i = 0; i < row_num_2; i++)
        bin_b[i] = new unsigned int[bin_col_num];
    //初始化为0
    for (int i = 0; i < row_num_1; i++)
    {
        for (int j = 0; j < bin_col_num; j++)
            bin_a[i][j] = 0;
    }
    for (int i = 0; i < row_num_2; i++)
    {
        for (int j = 0; j < bin_col_num; j++)
            bin_b[i][j] = 0;
    }
}
void dec_to_bin(int** x, unsigned int** bin_x, int row_num)
{
    for (int i = 0; i < row_num; i++)
    {
        int temp1;
        int temp2;
        int j = 0;
        while (x[i][j] != -1)
        {
            //bin_x[i][temp1]+=1(后面有temp2个0)->设置临时变量为1，并将其左移temp位
            temp1 = bin_col_num - 1 - x[i][j] / 32;
            temp2 = x[i][j] % 32;
            unsigned int temp3 = 1;
            temp3 <<= temp2;
            bin_x[i][temp1] += temp3;
            j++;
        }
    }
}
bool cmp(int x, int y)
{
    return x > y;
}
void bin_to_dec(int** x, vector<unsigned int*> bin_v)
{
    for (int i = 0; i < bin_v.size(); i++)
    {
        int num = 0;//每一行存储的数的个数
        for (int j = 0; j < bin_col_num; j++)
        {
            int mov_num = 0;//每个bin_x数组中的数向右移动的位数
            while (bin_v[i][j] != 0)
            {
                if (bin_v[i][j] & 1)//当前位的值为1
                {
                    int t = 32 * (bin_col_num - j - 1) + mov_num;//当前位置对应的十进制数的值
                    x[i][num] = t;
                    num++;//跳出循环后，该行存储了num个值
                }
                //右移
                bin_v[i][j] >>= 1;
                mov_num++;
            }
        }
        //对于这一行进行排序
        sort(x[i], x[i] + num, cmp);
        //将每行后加一个等于-1的位，作为结束的标识位
        x[i][num] = -1;
    }
}
void malloc_a1b1()
{
    a1 = new int* [bin_a_v.size()];
    b1 = new int* [bin_b_v.size()];
    for (int i = 0; i < bin_a_v.size(); i++)
        a1[i] = new int[col_num + 1];
    for (int i = 0; i < bin_b_v.size(); i++)
        b1[i] = new int[col_num + 1];
}

void d_mem()
{
    int nnn1 = bin_a_v.size();
    int nnn2 = bin_b_v.size();
    for (int i = 0; i < nnn1; i++)
    {
        delete[]a1[i];
    }
    for (int i = 0; i < nnn2; i++)
    {
        delete[]b1[i];
    }
    for (int i = 0; i < row_num_1; i++)
    {
        delete[]a[i];
        delete[]bin_a[i];
        //delete[]a1[i];
    }
    delete[]a;
    delete[]bin_a;
    delete[]a1;
    for (int i = 0; i < row_num_2; i++)
    {
        delete[]b[i];
        delete[]bin_b[i];
        //delete[]b1[i];
    }
    delete[]b;
    delete[]bin_b;
    delete[]b1;
}

int find_first(unsigned int* k)
{
    int pos_n = -1;
    for (int i = 0; i < bin_col_num; i++)
    {
        if (k[i] != 0)
        {
            unsigned int temp = k[i];
            int n = 0;
            //将该行右移，直到为0，则能确定第一个1所在位置
            while (temp != 0)
            {
                temp >>= 1;
                n++;
            }
            //pos_n = 32 * (i + 1) - n - 30;
            //cout <<"i:"<<i<<" " << "n:" << n << " " << "pos:" << pos_n << endl;
            pos_n = 32 * (bin_col_num - i - 1) + n - 1;
            break;
        }
    }
    return pos_n;
}

void handle_one(unsigned int* k1, unsigned int* k2, int n_bei, bool& judge1, bool& judge2)//k1->消元子,k2->被消元行
{
    /*
    //通过异或运算进行消元
    for (int i = 0; i < bin_col_num; i++)
       k2[i] = k1[i] ^ k2[i];
    */
    //可能的并行化操作
    int j;
    for(j=0;j+4<bin_col_num;j+=4)
    {
        //将k1,k2从内存加载到向量寄存器
        uint32x4_t v_k1=vld1q_u32(&k1[j]);
        uint32x4_t v_k2=vld1q_u32(&k2[j]);
        //向量对位异或
        v_k2=veorq_u32(v_k1,v_k2);
        //将数据从向量寄存器存储到内存
        vst1q_u32(&k2[j],v_k2);
    }
    for(;j<bin_col_num;j++)
        k2[j]=k1[j]^k2[j];
    
    //cout << find_first(k2) << endl;
    //消元后的处理
    //消元后：1）被消元行变为0->丢弃
    //        2）被消元行首项没有对应的消元子->升级为消元子（不再为被消元行）

    bool j1 = false;
    for (int i = 0; i < bin_col_num; i++)
    {
        if (k2[i] != 0)
        {
            j1 = true;
            break;
        }
    }
    if (j1 == false)//该被消元行为0->丢弃
    {
        //bin_b_v.erase(bin_b_v.begin() + n_bei, bin_b_v.begin() + n_bei + 1);
        //row_n_2--;
        judge1 = true;
        //sort(bin_b_v.begin(), bin_b_v.end(), cmp1);
    }
    int bei_first_1_loc = find_first(bin_b_v[n_bei]);

    if (equal_loc(bei_first_1_loc) == false)//被消元行首项没有对应的消元子
    {
        //把该被消元行bin_b_v[n_bei]变为消元子
        //找到该行应该插入到的位置
        int l = find_insert_loc(bei_first_1_loc);//消元子是按照首个非零元素所在位置从大到小排好的

        //bin_a_v.insert(bin_a_v.begin() + l, bin_b_v[n_bei]);
        unsigned int* temp = new unsigned[bin_col_num];
        for (int i = 0; i < bin_col_num; i++)
        {
            temp[i] = bin_b_v[n_bei][i];
        }
        bin_a_v.insert(bin_a_v.begin() + l, temp);
        //bin_a_v.insert(bin_a_v.begin()+l, bin_b_v[n_bei], bin_b_v[n_bei+1]);

        row_n_1++;
        judge2 = true;
    }
}
bool equal_loc(int first_1_loc)
{
    int judge = false;
    //遍历消元子
    for (int ui = 0; ui < bin_a_v.size(); ui++)
    {
        if (first_1_loc == find_first(bin_a_v[ui]))
        {
            judge = true;
            break;
        }
    }
    return judge;
}
int find_insert_loc(int& loc)
{
    int lll = -1;
    for (int i = 0; i < row_n_1; i++)
    {
        if (find_first(bin_a_v[i]) < loc)
        {
            lll = i;
            break;
        }
    }
    if (lll == -1)
        lll = row_n_1;
    return lll;
}