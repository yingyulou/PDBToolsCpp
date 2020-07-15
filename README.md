# PDBToolsCpp

PDB文件解析与坐标线性代数运算工具集。

PDB文件在PDBToolsCpp中将被解析为4个层级：Protein -> Chain -> Residue -> Atom

**本文档中所有的“角度”均指弧度制角度；所有的“旋转矩阵”均指右乘旋转矩阵。**

**任何结构对象指针均定义于堆内存上，使用时需注意内存管理。（详见下文析构函数部分内容）**

## 编译及使用说明

- 导入"PDBTools"头文件即可使用：

``` Cpp
#include <PDBToolsCpp/PDBTools>
```

- 编译依赖项：

1. boost
2. Eigen

- 链接依赖项：

1. boost_filesystem

- 编译器需支持C++11或以上标准

- PDBToolsCpp的所有接口均位于namespace PDBTools下

## PDB文件解析函数

### 1. Load

``` Cpp
Protein *Load(const string &pdbFilePath, bool parseHBool = false);
```

将PDB文件解析为Protein对象指针。

#### 参数：

- pdbFilePath：PDB文件路径
- parseHBool：是否开启氢原子解析

#### 返回值：

- Protein对象指针

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
```

### 2. LoadModel

``` Cpp
vector<Protein *> LoadModel(const string &pdbFilePath, bool parseHBool = false);
```

将含有"MODEL"关键词的PDB文件解析为Protein对象指针vector。

#### 参数：

- pdbFilePath：PDB文件路径
- parseHBool：是否开启氢原子解析

#### 返回值：

- Protein对象指针构成的vector

#### 例：

``` Cpp
vector<Protein *> proPtrList = LoadModel("xxxx.pdb");
```

## 属性

### Protein属性

#### 1. name

``` Cpp
string name;
```

PDB文件名（不包含".pdb"）。

#### 2. sub

``` Cpp
vector<Chain *> sub;
```

this包含的所有链对象指针。

### Chain属性

#### 1. name

``` Cpp
string name;
```

链名。

#### 2. owner

``` Cpp
Protein *owner;
```

this所属的Protein。

#### 3. sub

``` Cpp
vector<Residue *> sub;
```

this包含的所有残基对象指针。

### Residue属性

#### 1. name

``` Cpp
string name;
```

残基名。

#### 2. num

``` Cpp
int num;
```

残基编号。

#### 3. ins

``` Cpp
string ins;
```

残基插入字符。

#### 4. owner

``` Cpp
Chain *owner;
```

this所属的Chain。

#### 5. sub

``` Cpp
vector<Atom *> sub;
```

this包含的所有原子对象指针。

### Atom属性

#### 1. name

``` Cpp
string name;
```

原子名。

#### 2. num

``` Cpp
int num;
```

原子编号。

#### 3. coord

``` Cpp
Eigen::RowVector3d coord;
```

原子坐标。

#### 4. alt

``` Cpp
string alt;
```

备用位置指示符。

#### 5. occ

``` Cpp
string occ;
```

占有。

#### 6. tempF

``` Cpp
string tempF;
```

温度因子。

#### 7. ele

``` Cpp
string ele;
```

元素符号。

#### 8. chg

``` Cpp
string chg;
```

电荷。

#### 9. owner

``` Cpp
Residue *owner;
```

this所属的Residue。

## 成员函数

### 构造函数

#### 1. Protein构造函数

``` Cpp
explicit Protein(const string &proteinID = "");
```

#### 参数：

- proteinID：蛋白名，用于初始化name属性

#### 例：

``` Cpp
Protein *proPtr = new Protein("xxxx");
```

#### 2. Chain构造函数

``` Cpp
explicit Chain(const string &chainName = "", Protein *chainOwner = nullptr);
```

#### 参数：

- chainName：链名，用于初始化name属性
- chainOwner：this的所属Protein，用于初始化owner属性。如果owner不为nullptr，则构造函数将自动在owner与this之间建立从属关系

#### 例：

``` Cpp
Protein *proPtr = new Protein("xxxx");
Chain *chainPtr = new Chain("X", proPtr);
```

#### 3. Residue构造函数

``` Cpp
explicit Residue(const string &resName = "", int resNum = 0,
    const string &resIns = "", Chain *resOwner = nullptr);
```

#### 参数：

- resName：残基名，用于初始化name属性
- resNum：残基编号，用于初始化num属性
- resIns：残基插入字符，用于初始化ins属性
- resOwner：this的所属Chain，用于初始化owner属性。如果owner不为nullptr，则构造函数将自动在owner与this之间建立从属关系

#### 例：

``` Cpp
Chain *chainPtr = new Chain("X");
Residue *resPtr = new Residue("XXX", 0, "", chainPtr);
```

#### 4. Atom构造函数

``` Cpp
explicit Atom(const string &atomName = "", int atomNum = 0,
    const Eigen::RowVector3d &atomCoord = Eigen::RowVector3d(0., 0., 0.),
    const string &atomAltLoc = "", const string &atomOccupancy = "",
    const string &atomTempFactor = "", const string &atomElement = "",
    const string &atomCharge = "", Residue *atomOwner = nullptr);
```

#### 参数：

- atomName：原子名，用于初始化name属性
- atomNum：原子编号，用于初始化num属性
- atomCoord：原子坐标，用于初始化coord属性
- atomAltLoc：备用位置指示符，用于初始化alt属性
- atomOccupancy：占有，用于初始化occ属性
- atomTempFactor：温度因子，用于初始化tempF属性
- atomElement：元素符号，用于初始化ele属性
- atomCharge：电荷，用于初始化chg属性
- atomOwner：this的所属Residue，用于初始化owner属性。如果owner不为nullptr，则构造函数将自动在owner与this之间建立从属关系

#### 例：

``` Cpp
Residue *resPtr = new Residue("X");
Atom *atomPtr = new Atom("X", 0, RowVector3d(0., 0., 0.), "", "", "", "", "", resPtr);
```

### 析构函数

任何结构对象的析构函数将递归地析构此结构对象所属的一切子结构对象。

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
delete proPtr;
```

### 迭代器

非Atom层级均支持委托至sub属性的迭代器：

``` Cpp
typedef typename vector<SubType *>::iterator iterator;
typedef typename vector<SubType *>::const_iterator const_iterator;
typedef typename vector<SubType *>::reverse_iterator reverse_iterator;
typedef typename vector<SubType *>::const_reverse_iterator const_reverse_iterator;

iterator begin();
const_iterator begin() const;
iterator end();
const_iterator end() const;
const_iterator cbegin() const;
const_iterator cend() const;
reverse_iterator rbegin();
const_reverse_iterator rbegin() const;
reverse_iterator rend();
const_reverse_iterator rend() const;
const_reverse_iterator crbegin() const;
const_reverse_iterator crend() const;
```

#### 例：

``` Cpp
Protein *proPtr = Load("6urw.pdb");

for (Chain *chainPtr: *proPtr)
{
    for (Residue *resPtr: *chainPtr)
    {
        for (Atom *atomPtr: *resPtr)
        {
            cout << *atomPtr << endl;
        }
    }
}
```

### 所有层级公有成员函数

#### 1. str

``` Cpp
string str() const;
```

得到结构对象的摘要信息。

#### 参数：

- void

#### 返回值：

- 结构对象的摘要信息

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
cout << proPtr->str() << endl;  // 相当于cout << *proPtr << endl;
```

#### 2. Dump

``` Cpp
SelfType *Dump(const string &dumpFilePath, const string &fileMode = "w");

const SelfType *Dump(const string &dumpFilePath, const string &fileMode = "w") const;
```

将this输出到PDB文件。

#### 参数：

- dumpFilePath：输出PDB文件路径
- fileMode：文件句柄打开模式

#### 返回值：

- this

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
proPtr->Dump("xxxx.pdb");
```

#### 3. Dumps

``` Cpp
string Dumps() const;
```

得到字符串形式的PDB文件内容。

#### 参数：

- void

#### 返回值：

- 字符串形式的PDB文件内容

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
string dumpStr = proPtr->Dumps();
```

#### 4. Copy

``` Cpp
Protein *Copy() const;
Chain *Copy() const;
Residue *Copy() const;
Atom *Copy() const;
```

得到this的深拷贝。

#### 参数：

- void

#### 返回值：

- this的深拷贝

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Protein *copyProPtr = proPtr->Copy();
```

### 非Atom层级公有成员函数

#### 1. GetResidues

``` Cpp
vector<Residue *> GetResidues();

vector<const Residue *> GetResidues() const;
```

跨层级直接返回this包含的所有残基对象指针。

#### 参数：

- void

#### 返回值：

- this包含的所有残基对象指针

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
vector<Residue *> resPtrList = proPtr->GetResidues();
```

#### 2. GetAtoms, FilterAtoms, GetAtomsCoord, FilterAtomsCoord

``` Cpp
vector<Atom *> GetAtoms();

vector<const Atom *> GetAtoms() const;

vector<Atom *> FilterAtoms(const unordered_set<string> &atomNameSet = {"CA"});

vector<const Atom *> FilterAtoms(const unordered_set<string> &atomNameSet = {"CA"}) const;

vector<Eigen::RowVector3d *> GetAtomsCoord();

vector<const Eigen::RowVector3d *> GetAtomsCoord() const;

vector<Eigen::RowVector3d *> FilterAtomsCoord(
    const unordered_set<string> &atomNameSet = {"CA"});

vector<const Eigen::RowVector3d *> FilterAtomsCoord(
    const unordered_set<string> &atomNameSet = {"CA"}) const;
```

跨层级直接返回this包含的所有，或按原子的name属性筛选后的原子对象指针或原子坐标指针。

#### 参数：

- atomNameSet：原子名集合

#### 返回值：

- 原子对象指针列表/原子坐标指针列表

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
vector<Atom *> atomPtrList = proPtr->GetAtoms();
vector<Atom *> filterAtomPtrList = proPtr->FilterAtoms({"N", "CA", "C"});
vector<RowVector3d *> atomCoordList = proPtr->GetAtomsCoord();
vector<RowVector3d *> filterAtomCoordList = proPtr->FilterAtomsCoord({"N", "CA", "C"});
```

#### 3. center

``` Cpp
Eigen::RowVector3d center() const;
```

得到this的所有原子坐标的几何中心。

#### 参数：

- void

#### 返回值：

- this的所有原子坐标的几何中心

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
RowVector3d centerCoord = proPtr->center();
```

#### 4. MoveCenter

``` Cpp
SelfType *MoveCenter();
```

将this的所有原子坐标减去center()向量。

#### 参数：

- void

#### 返回值：

- this

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
proPtr->MoveCenter();
```

#### 5. seq

``` Cpp
string seq() const;
```

得到this的残基序列。

#### 参数：

- void

#### 返回值：

- this的残基序列

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
string seqStr = proPtr->seq();
```

#### 6. fasta

``` Cpp
string fasta(const string &titleStr = "") const;
```

得到字符串形式的Fasta文件内容。

#### 参数：

- titleStr：Fasta文件的标题内容，如果传入空字符串，则标题将被设置为this->name

#### 返回值：

- 字符串形式的Fasta文件内容

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
string fastaStr = proPtr->fasta("xxxx");
```

#### 7. DumpFasta

``` Cpp
SelfType *DumpFasta(const string &dumpFilePath,
    const string &titleStr = "", const string &fileMode = "w");

const SelfType *DumpFasta(const string &dumpFilePath,
    const string &titleStr = "", const string &fileMode = "w") const;
```

将this输出到Fasta文件。

#### 参数：

- dumpFilePath：输出Fasta文件路径
- titleStr：Fasta文件的标题内容，如果传入空字符串，则标题将被设置为this->name
- fileMode：文件句柄打开模式

#### 返回值：

- this

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
proPtr->DumpFasta("xxxx.fasta", "xxxx");
```

#### 8. RenumResidues, RenumAtoms

``` Cpp
SelfType *RenumResidues(int startNum = 1);

SelfType *RenumAtoms(int startNum = 1);
```

对this的所有残基/原子进行重编号。

#### 参数：

- startNum：起始编号

#### 返回值：

- this

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
proPtr->RenumResidues()->RenumAtoms();
```

#### 9. PushBack, Insert

``` Cpp
SelfType *PushBack(const SubType *subPtr);

SelfType *Insert(const_iterator insertIter, SubType *subPtr);
SelfType *Insert(const_iterator insertIter, int n, SubType *subPtr);

template <typename ForwardIterator>
SelfType *Insert(const_iterator insertIter, ForwardIterator firstIter,
    ForwardIterator lastIter);

SelfType *Insert(const_iterator insertIter, initializer_list<SubType *> initializerList);
```

为this追加/插入子结构。

所有添加至this的子结构都是原结构对象指针调用Copy成员函数得到的拷贝，且会与this建立从属关系。

#### 参数：

- subPtr：this对应的子结构对象指针
- insertIter：插入位置迭代器
- n：subPtr的重复次数
- firstIter：this对应的子结构左迭代器
- lastIter：this对应的子结构右迭代器
- initializerList：this对应的子结构对象指针列表

#### 返回值：

- this

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");

proPtr->
    PushBack(proPtr->sub[0])->
    Insert(proPtr->begin(), proPtr->sub[0])->
    Insert(proPtr->begin(), 2, proPtr->sub[0])->
    Insert(proPtr->begin(), proPtr->begin(), proPtr->end())->
    Insert(proPtr->begin(), {proPtr->sub[0], proPtr->sub[0]});
```

#### 10. MoveBack, MoveInsert

``` Cpp
SelfType *MoveBack(SubType *subPtr);

SelfType *MoveInsert(const_iterator insertIter, SubType *subPtr);
```

为this移动并追加/移动并插入子结构。

移动至this的子结构会断开其原从属关系（如果有），并与this建立新的从属关系。

#### 参数：

- subPtr：this对应的子结构对象指针
- insertIter：插入位置迭代器

#### 返回值：

- this

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");

proPtr->
    MoveBack(proPtr->sub[0])->
    MoveInsert(proPtr->begin(), proPtr->sub[0]);
```

#### 11. RemoveAlt

``` Cpp
SelfType *RemoveAlt();
```

遍历this包含的所有原子对象指针，如果原子对象指针的alt属性为""，则忽略，如果为"A"，则修改为""，否则删除当前原子。

#### 参数：

- void

#### 返回值：

- this

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
proPtr->RemoveAlt();
```

#### 12. subMap

``` Cpp
// Protein
unordered_map<string, Chain *> subMap();

unordered_map<string, const Chain *> subMap() const;


// Chain
unordered_map<string, Residue *> subMap();

unordered_map<string, const Residue *> subMap() const;


// Residue
unordered_map<string, Atom *> subMap();

unordered_map<string, const Atom *> subMap() const;
```

对于Protein对象：得到this包含的所有链名 -> 链对象指针哈希表。

对于Chain对象：得到this包含的所有完整残基名 -> 残基对象指针哈希表。

对于Residue对象：得到this包含的所有原子名 -> 原子对象指针哈希表。

#### 参数：

- void

#### 返回值：

- 对于Protein对象：this包含的所有链名 -> 链对象指针哈希表
- 对于Chain对象：this包含的所有完整残基名 -> 残基对象指针哈希表
- 对于Residue对象：this包含的所有原子名 -> 原子对象指针哈希表

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
unordered_map<string, Chain *> proSubMap = proPtr->subMap();
unordered_map<string, Residue *> chainSubMap = proPtr->sub[0]->subMap();
unordered_map<string, Atom *> resSubMap = proPtr->sub[0]->sub[0]->subMap();
```

### 非Protein层级公有成员函数

#### 1. iter

``` Cpp
typename vector<SelfType *>::iterator iter();
typename vector<SelfType *>::const_iterator iter() const;
```

得到this在this->owner->sub中的位置迭代器。

#### 参数：

- void

#### 返回值：

- this在this->owner->sub中的位置迭代器

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
vector<Chain *>::iterator iterInOwner = proPtr->sub[0]->iter();  // begin
```

#### 2. pre, next

``` Cpp
SelfType *pre(int shiftLen = 1);

const SelfType *pre(int shiftLen = 1) const;

SelfType *next(int shiftLen = 1);

const SelfType *next(int shiftLen = 1) const;
```

得到this在this->owner->sub中的前/后第N个同级结构对象指针

#### 参数：

- shiftLen：偏移量

#### 返回值：

- this的前/后第N个同级结构对象指针

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Chain *chainPtr = proPtr->sub[1];
Chain *preChainPtr = chainPtr->pre();
Chain *nextChainPtr = chainPtr->next();
chainPtr->pre(2);  // Error!
```

#### 3. Remove

``` Cpp
typename vector<SelfType *>::iterator Remove(bool deteleBool = true);
```

从this->owner->sub中删除并析构this。

#### 参数：

- deteleBool：如果deteleBool参数为false，则不会析构this

#### 返回值：

- 删除this后，this的后继元素位置迭代器

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
proPtr->sub[0]->Remove();
```

### Protein的其他成员函数

无。

### Chain的其他成员函数

无。

### Residue的其他成员函数

#### 1. compNum

``` Cpp
string compNum() const;

Residue *compNum(int resNum, const string &resIns = "");
```

同时获取/设定残基对象的num + ins属性。

#### 参数：

- resNum：残基编号
- resIns：残基插入字符

#### 返回值：

- 残基对象的num + ins属性字符串（对于无参版本）/this（对于有参版本）

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Residue *resPtr = proPtr->sub[0]->sub[0];
resPtr->compNum(0, "");
string compNum = resPtr->compNum();
```

#### 2. coordMap

``` Cpp
unordered_map<string, Eigen::RowVector3d *> coordMap();

unordered_map<string, const Eigen::RowVector3d *> coordMap() const;
```

得到this包含的所有原子名 -> 原子坐标指针哈希表。

#### 参数：

- void

#### 返回值：

- this包含的所有原子名 -> 原子坐标指针哈希表

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Residue *resPtr = proPtr->sub[0]->sub[0];
unordered_map<string, RowVector3d *> coordMap = resPtr->coordMap();
```

#### 3. 二面角相关

详见下文。

### Atom的其他成员函数

#### 1. operator-

``` Cpp
double operator-(const Atom &rhs) const;
```

计算两原子间的欧几里得距离。

#### 参数：

- rhs：另一个原子对象

#### 返回值：

- 两原子间的欧几里得距离

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
double atomDis = *proPtr->sub[0]->sub[0]->sub[0] - *proPtr->sub[0]->sub[0]->sub[1];
```

## 残基二面角

残基对象实现了若干对蛋白主/侧链二面角进行计算和旋转相关的成员函数（即以下所有成员函数的this都专指Residue *）。

### 主链二面角

**对主链进行操作时请注意：N端与C端的两个残基分别无法进行二面角Phi与Psi的计算或调整（因为这两个二面角不存在）。如果出现上述情况，则将抛出out_of_range异常。**

#### 1. CalcBBDihedralAngle

``` Cpp
double CalcBBDihedralAngle(DIH dihedralEnum) const;
```

计算主链二面角。

#### 参数：

- dihedralEnum：主链二面角种类。DIH::PHI或DIH::L表示Phi；DIH::PSI或DIH::R表示Psi

#### 返回值：

- 有符号二面角值（-pi ~ pi）

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Residue *resPtr = proPtr->sub[0]->sub[1];
double dihedralAngle = resPtr->CalcBBDihedralAngle(DIH::PHI);
```

#### 2. CalcBBRotationMatrixByDeltaAngle, CalcBBRotationMatrixByTargetAngle

``` Cpp
Residue *CalcBBRotationMatrixByDeltaAngle(DIH dihedralEnum,
    SIDE sideEnum, double deltaAngle, Eigen::RowVector3d &moveCoord,
    Eigen::Matrix3d &rotationMatrix);

const Residue *CalcBBRotationMatrixByDeltaAngle(DIH dihedralEnum,
    SIDE sideEnum, double deltaAngle, Eigen::RowVector3d &moveCoord,
    Eigen::Matrix3d &rotationMatrix) const;

Residue *CalcBBRotationMatrixByTargetAngle(DIH dihedralEnum,
    SIDE sideEnum, double targetAngle, Eigen::RowVector3d &moveCoord,
    Eigen::Matrix3d &rotationMatrix);

const Residue *CalcBBRotationMatrixByTargetAngle(DIH dihedralEnum,
    SIDE sideEnum, double targetAngle, Eigen::RowVector3d &moveCoord,
    Eigen::Matrix3d &rotationMatrix) const;
```

以旋转角度/目标角度作为参数，计算主链旋转矩阵。

#### 参数：

- dihedralEnum：主链二面角种类。DIH::PHI或DIH::L表示Phi；DIH::PSI或DIH::R表示Psi
- sideEnum：转动侧。SIDE::N或SIDE::L表示转动N端；SIDE::C或SIDE::R表示转动C端
- deltaAngle/targetAngle：旋转角度/目标角度
- moveCoord：旋转前/后平移向量（返回值）
- rotationMatrix：旋转矩阵（返回值）

#### 返回值：

- this

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Residue *resPtr = proPtr->sub[0]->sub[1];
RowVector3d moveCoord;
Matrix3d rotationMatrix;
resPtr->CalcBBRotationMatrixByDeltaAngle(DIH::PHI, SIDE::N, 1., moveCoord, rotationMatrix);
resPtr->CalcBBRotationMatrixByTargetAngle(DIH::PHI, SIDE::N, 0., moveCoord, rotationMatrix);
```

#### 3. GetBBRotationAtomObj

``` Cpp
vector<Atom *> GetBBRotationAtomPtr(DIH dihedralEnum, SIDE sideEnum);

vector<const Atom *> GetBBRotationAtomPtr(DIH dihedralEnum, SIDE sideEnum) const;
```

获取以给定参数进行旋转时，所有需要旋转的原子对象指针列表。

#### 参数：

- dihedralEnum：主链二面角种类。DIH::PHI或DIH::L表示Phi；DIH::PSI或DIH::R表示Psi
- sideEnum：转动侧。SIDE::N或SIDE::L表示转动N端；SIDE::C或SIDE::R表示转动C端

#### 返回值：

- 所有需要旋转的原子对象指针列表

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Residue *resPtr = proPtr->sub[0]->sub[1];
vector<Atom *> rotationAtomPtrList = resPtr->GetBBRotationAtomPtr(DIH::PHI, SIDE::N);
```

#### 4. RotateBBDihedralAngleByDeltaAngle, RotateBBDihedralAngleByTargetAngle

``` Cpp
Residue *RotateBBDihedralAngleByDeltaAngle(DIH dihedralEnum,
    SIDE sideEnum, double deltaAngle);

Residue *RotateBBDihedralAngleByTargetAngle(DIH dihedralEnum,
    SIDE sideEnum, double targetAngle);
```

以旋转角度/目标角度作为参数直接旋转主链。

#### 参数：

- dihedralEnum：主链二面角种类。DIH::PHI或DIH::L表示Phi；DIH::PSI或DIH::R表示Psi
- sideEnum：转动侧。SIDE::N或SIDE::L表示转动N端；SIDE::C或SIDE::R表示转动C端
- deltaAngle/targetAngle：旋转角度/目标角度

#### 返回值：

- this

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Residue *resPtr = proPtr->sub[0]->sub[1];
resPtr->RotateBBDihedralAngleByDeltaAngle(DIH::PHI, SIDE::N, 1.);
resPtr->RotateBBDihedralAngleByTargetAngle(DIH::PHI, SIDE::N, 0.);
```

### 侧链二面角

**对侧链进行调整时请注意：GLY、ALA残基由于不存在侧链二面角，不可调用下列成员函数。且不可使用不存在的侧链二面角索引值调用下列成员函数。如果出现上述情况，则将抛出out_of_range异常。**

#### 1. CalcSCDihedralAngle

``` Cpp
double CalcSCDihedralAngle(int dihedralIdx) const;
```

计算侧链二面角。

#### 参数：

- dihedralIdx：侧链二面角索引值。索引值从0开始编号，最大允许索引值根据残基种类而不同。索引值表示某个残基从主链到侧链方向上的第N个侧链二面角

#### 返回值：

- 有符号二面角值（-pi ~ pi）

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Residue *resPtr = proPtr->sub[0]->sub[1];
double dihedralAngle = resPtr->CalcSCDihedralAngle(0);
```

#### 2. CalcSCRotationMatrixByDeltaAngle, CalcSCRotationMatrixByDeltaAngle

``` Cpp
Residue *CalcSCRotationMatrixByDeltaAngle(int dihedralIdx, double deltaAngle,
    Eigen::RowVector3d &moveCoord, Eigen::Matrix3d &rotationMatrix);

const Residue *CalcSCRotationMatrixByDeltaAngle(int dihedralIdx, double deltaAngle,
    Eigen::RowVector3d &moveCoord, Eigen::Matrix3d &rotationMatrix) const;

Residue *CalcSCRotationMatrixByTargetAngle(int dihedralIdx, double targetAngle,
    Eigen::RowVector3d &moveCoord, Eigen::Matrix3d &rotationMatrix);

const Residue *CalcSCRotationMatrixByTargetAngle(int dihedralIdx, double targetAngle,
    Eigen::RowVector3d &moveCoord, Eigen::Matrix3d &rotationMatrix) const;
```

以旋转角度/目标角度作为参数，计算侧链旋转矩阵。

#### 参数：

- dihedralIdx：侧链二面角索引值。索引值从0开始编号，最大允许索引值根据残基种类而不同。索引值表示某个残基从主链到侧链方向上的第N个侧链二面角
- deltaAngle/targetAngle：旋转角度/目标角度
- moveCoord：旋转前/后平移向量（返回值）
- rotationMatrix：旋转矩阵（返回值）

#### 返回值：

- this

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Residue *resPtr = proPtr->sub[0]->sub[1];
RowVector3d moveCoord;
Matrix3d rotationMatrix;
resPtr->CalcSCRotationMatrixByDeltaAngle(0, 1., moveCoord, rotationMatrix);
resPtr->CalcSCRotationMatrixByTargetAngle(0, 0., moveCoord, rotationMatrix);
```

#### 3. GetSCRotationAtomObj

``` Cpp
vector<Atom *> GetSCRotationAtomPtr(int dihedralIdx);

vector<const Atom *> GetSCRotationAtomPtr(int dihedralIdx) const;
```

获取以给定侧链二面角进行旋转时，所有需要旋转的原子对象指针列表。

#### 参数：

- dihedralIdx：侧链二面角索引值。索引值从0开始编号，最大允许索引值根据残基种类而不同。索引值表示某个残基从主链到侧链方向上的第N个侧链二面角

#### 返回值：

- 所有需要旋转的原子对象指针列表

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Residue *resPtr = proPtr->sub[0]->sub[1];
vector<Atom *> rotationAtomPtrList = resPtr->GetSCRotationAtomPtr(0);
```

#### 4. RotateSCDihedralAngleByDeltaAngle, RotateSCDihedralAngleByTargetAngle

``` Cpp
Residue *RotateSCDihedralAngleByDeltaAngle(int dihedralIdx, double deltaAngle);

Residue *RotateSCDihedralAngleByTargetAngle(int dihedralIdx, double targetAngle);
```

以旋转角度/目标角度作为参数直接旋转侧链。

#### 参数：

- dihedralIdx：侧链二面角索引值。索引值从0开始编号，最大允许索引值根据残基种类而不同。索引值表示某个残基从主链到侧链方向上的第N个侧链二面角
- deltaAngle/targetAngle：旋转角度/目标角度

#### 返回值：

- this

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
Residue *resPtr = proPtr->sub[0]->sub[1];
resPtr->RotateSCDihedralAngleByDeltaAngle(0, 1.);
resPtr->RotateSCDihedralAngleByTargetAngle(0, 0.);
```

## 数学函数

### 1. Degrees, Radians

``` Cpp
double Degrees(double radiansAngle);
double Radians(double degreesAngle);
```

角度，弧度互相转换。

#### 参数：

- radiansAngle，degreesAngle：弧度值/角度值

#### 返回值：

- 转换后的角度值/弧度值

#### 例：

``` Cpp
cout << Degrees(1.) << endl;
cout << Radians(1.) << endl;
```

### 2. CalcVectorAngle

``` Cpp
double CalcVectorAngle(const Eigen::RowVector3d &coordA, const Eigen::RowVector3d &coordB);
```

计算两向量夹角。

#### 参数：

- coordA，coordB：两个三维向量

#### 返回值：

- 两向量夹角值（0 ~ pi）

#### 例：

``` Cpp
double vectorAngle = CalcVectorAngle(RowVector3d(1., 2., 3.), RowVector3d(4., 5., 6.));
```

### 3. CalcRotationMatrix

``` Cpp
Eigen::Matrix3d CalcRotationMatrix(const Eigen::RowVector3d &rotationAxis, double rotationAngle);
```

计算轴角旋转矩阵。

#### 参数：

- rotationAxis：旋转轴向量，无需缩放至单位长度
- rotationAngle：旋转角

#### 返回值：

- 旋转矩阵

#### 例：

``` Cpp
Matrix3d rotationMatrix = CalcRotationMatrix(RowVector3d(1., 2., 3.), 1.);
```

### 4. CalcRotationMatrixByTwoVector

``` Cpp
Eigen::Matrix3d CalcRotationMatrixByTwoVector(
    const Eigen::RowVector3d &coordA, const Eigen::RowVector3d &coordB);
```

计算从向量A旋转至向量B所需要的旋转矩阵。

#### 参数：

- coordA，coordB：两个三维向量

#### 返回值：

- 旋转矩阵

#### 例：

``` Cpp
Matrix3d rotationMatrix = CalcRotationMatrixByTwoVector(
    RowVector3d(1., 2., 3.), RowVector3d(4., 5., 6.));
```

### 5. CalcDihedralAngle

``` Cpp
double CalcDihedralAngle(
    const Eigen::RowVector3d &coordA, const Eigen::RowVector3d &coordB,
    const Eigen::RowVector3d &coordC, const Eigen::RowVector3d &coordD);
```

计算二面角。

#### 参数：

- coordA，coordB，coordC，coordD：四个三维向量

#### 返回值：

- 有符号二面角值（-pi ~ pi）

#### 例：

``` Cpp
double dihedralAngle = CalcDihedralAngle(RowVector3d(1., 2., 3.), RowVector3d(4., 5., 6.),
    RowVector3d(7., 8., 9.), RowVector3d(10., 11., 12.));
```

### 6. CalcRMSD

``` Cpp
double CalcRMSD(const Eigen::Matrix<double, Eigen::Dynamic, 3> &coordArrayA,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &coordArrayB);
```

对两组等长的三维坐标计算RMSD。

#### 参数：

- coordArrayA，coordArrayB：两组等长的矩阵（N * 3）

#### 返回值：

- RMSD值

#### 例：

``` Cpp
Matrix<double, 2, 3> coordArrayA, coordArrayB;
coordArrayA << 1., 2., 3., 4., 5., 6.;
coordArrayB << 7., 8., 9., 10., 11., 12.;
double rmsdValue = CalcRMSD(coordArrayA, coordArrayB);
```

### 7. CalcSuperimposeRotationMatrix

``` Cpp
void CalcSuperimposeRotationMatrix(
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &sourceCoordArray,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &targetCoordArray,
    Eigen::RowVector3d &sourceCenterCoord,
    Eigen::Matrix3d &rotationMatrix,
    Eigen::RowVector3d &targetCenterCoord);
```

计算从sourceCoordArray到targetCoordArray的叠合旋转矩阵。

此函数将使得((sourceCoordArray.rowwise() - sourceCenterCoord) * rotationMatrix).rowwise() + targetCenterCoord与targetCoordArray形成叠合（RMSD最小）。

#### 参数：

- sourceCoordArray, targetCoordArray：两组等长的矩阵（N * 3）
- sourceCenterCoord, targetCenterCoord：前/后平移向量（返回值）
- rotationMatrix：旋转矩阵（返回值）

#### 返回值：

- void

#### 例：

``` Cpp
Matrix<double, 2, 3> coordArrayA, coordArrayB;
coordArrayA << 1., 2., 3., 4., 5., 6.;
coordArrayB << 7., 8., 9., 10., 11., 12.;

RowVector3d sourceCenterCoord, targetCenterCoord;
Matrix3d rotationMatrix;

CalcSuperimposeRotationMatrix(coordArrayA, coordArrayB,
    sourceCenterCoord, rotationMatrix, targetCenterCoord);

cout << ((coordArrayA.rowwise() - sourceCenterCoord) * rotationMatrix).rowwise() +
    targetCenterCoord << endl << coordArrayB << endl;
```

### 8. CalcRMSDAfterSuperimpose

``` Cpp
double CalcRMSDAfterSuperimpose(
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &coordArrayA,
    const Eigen::Matrix<double, Eigen::Dynamic, 3> &coordArrayB);
```

叠合并计算RMSD。

此函数会将coordArrayA通过CalcSuperimposeRotationMatrix函数向coordArrayB进行叠合，然后计算两组坐标之间的RMSD。

#### 参数：

- coordArrayA，coordArrayB：两组等长的矩阵（N * 3）

#### 返回值：

- RMSD值

#### 例：

``` Cpp
Matrix<double, 2, 3> coordArrayA, coordArrayB;
coordArrayA << 1., 2., 3., 4., 5., 6.;
coordArrayB << 7., 8., 9., 10., 11., 12.;
double rmsdValue = CalcRMSDAfterSuperimpose(coordArrayA, coordArrayB);
```

### 9. Coord2Matrix

``` Cpp
Eigen::Matrix<double, Eigen::Dynamic, 3> Coord2Matrix(
    const vector<Eigen::RowVector3d *> &coordPtrList);

Eigen::Matrix<double, Eigen::Dynamic, 3> Coord2Matrix(
    const vector<const Eigen::RowVector3d *> &coordPtrList);
```

将GetAtomsCoord或FilterAtomsCoord函数的返回值（vector\<Eigen::RowVector3d *\>类型）转为Eigen::Matrix\<double, Dynamic, 3\>类型。

#### 参数：

- coordPtrList：类型为vector\<RowVector3d *\>的vector

#### 返回值：

- 由coordPtrList参数转换得到的Matrix\<double, Dynamic, 3\>类型矩阵

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
vector<RowVector3d *> coordPtrList = proPtr->GetAtomsCoord();
Matrix<double, Dynamic, 3> coordMatrix = Coord2Matrix(coordPtrList);
```

## 常量

### 1. DIH

``` Cpp
enum class DIH;
```

枚举变量，表示主链二面角种类。DIH::PHI或DIH::L表示Phi，DIH::PSI或DIH::R表示Psi。

### 2. SIDE

``` Cpp
enum class SIDE;
```

枚举变量，表示主链二面角旋转时的转动侧。SIDE::N或SIDE::L表示转动N端，SIDE::C或SIDE::R表示转动C端。

### 3. RESIDUE_NAME_THREE_TO_ONE_MAP, RESIDUE_NAME_ONE_TO_THREE_MAP

``` Cpp
const unordered_map<string, string> RESIDUE_NAME_THREE_TO_ONE_MAP;
const unordered_map<string, string> RESIDUE_NAME_ONE_TO_THREE_MAP;
```

三字母，单字母残基名的相互转换哈希表。

## 其他函数

### 1. operator<<

``` Cpp
template <typename SelfType>
ostream &operator<<(ostream &os, const __StructBase<SelfType> &structObj);
```

任何结构对象均可用通过operator\<\<输出此对象的信息摘要。

#### 参数：

- os：输出流对象
- structPtr：任何层级对象

#### 返回值：

- 输出流对象

#### 例：

``` Cpp
Protein *proPtr = Load("xxxx.pdb");
cout << *proPtr << endl;
```

### 2. IsH

``` Cpp
inline bool IsH(const string &atomName);
```

判断一个原子名是否为氢原子。

#### 参数：

- atomName：原子名

#### 返回值：

- 原子名是否为氢原子

#### 例：

``` Cpp
bool isHBool = IsH("1H");
```

### 3. SplitCompNum

``` Cpp
inline void SplitCompNum(const string &compNumStr, int &resNum, string &resIns);
```

将完整残基编号分割为残基编号和残基插入编号。

#### 参数：

- compNumStr：完整残基编号
- resNum：残基编号（返回值）
- resIns：残基插入编号（返回值）

#### 返回值：

- void

#### 例：

``` Cpp
int resNum;
string resIns;
SplitCompNum("1A", resNum, resIns);
```

### 4. Dumpl

``` Cpp
template <typename SelfType>
void Dumpl(const vector<SelfType *> &structPtrList,
    const string &dumpFilePath, const string &fileMode = "w");
```

将任何对象指针构成的vector输出到PDB文件。

#### 参数：

- structPtrList：任何层级对象指针构成的vector
- dumpFilePath：输出PDB文件路径
- fileMode：文件句柄打开模式

#### 返回值：

- void

#### 例：

``` Cpp
vector<Protein *> proPtrList = LoadModel("xxxx.pdb");
Dumpl(proPtrList, "xxxx.pdb");
```

### 5. Dumpls

``` Cpp
template <typename SelfType>
string Dumpls(const vector<SelfType *> &structPtrList);
```

得到字符串形式的Dumpl函数输出内容。

#### 参数：

- structPtrList：任何层级对象指针构成的vector

#### 返回值：

- 字符串形式的Dumpl函数输出内容

#### 例：

``` Cpp
vector<Protein *> proPtrList = LoadModel("xxxx.pdb");
string dumpStr = Dumpls(proPtrList);
```

### 6. DumpFastal

``` Cpp
template <typename SelfType>
void DumpFastal(const vector<SelfType *> &structPtrList,
    const string &dumpFilePath, const string &fileMode = "w");
```

将非Atom对象指针构成的vector输出到Fasta文件。

#### 参数：

- structPtrList：非Atom对象指针构成的vector
- dumpFilePath：输出Fasta文件路径
- fileMode：文件句柄打开模式

#### 返回值：

- void

#### 例：

``` Cpp
vector<Protein *> proPtrList = LoadModel("xxxx.pdb");
DumpFastal(proPtrList, "xxxx.fasta");
```

### 7. DumpFastals

``` Cpp
template <typename SelfType>
string DumpFastals(const vector<SelfType *> &structPtrList);
```

得到字符串形式的DumpFastal函数输出内容。

#### 参数：

- structPtrList：非Atom对象指针构成的vector

#### 返回值：

- 字符串形式的DumpFastal函数输出内容

#### 例：

``` Cpp
vector<Protein *> proPtrList = LoadModel("xxxx.pdb");
string dumpStr = DumpFastals(proPtrList);
```

## 补充说明

### 解析函数

- PDB文件解析函数（Load、LoadModel）将完全按照PDB文件内容进行解析，不会对结构进行任何排序、合并或重组操作
- Load函数在解析时会跳过任何非"ATOM"关键词开头的行（包括"MODEL"）；而LoadModel函数会跳过任何非"ATOM"或"MODEL"关键词开头的行
- 解析时会去除所有字符串类型属性双端的空格字符

### 对于创建新对象指针的判定

#### 1. Load函数：

- Protein：只会在解析开始前创建唯一的一个，并最终返回这个对象指针
- Chain：解析开始时，以及每次检测到链名发生变化时（从上一个"ATOM"行到当前行），都会创建一个新的链对象指针
- Residue：解析开始时，创建新链时，以及每次检测到残基名、残基编号或残基插入字符三者之一发生变化时（从上一个"ATOM"行到当前行），都会创建一个新的残基对象指针
- Atom：每检测到一个新的"ATOM"行都会创建一个新的Atom对象指针

#### 2. LoadModel函数：

- Protein：解析开始前，以及每次检测到"MODEL"关键词时，都会创建一个新的蛋白对象指针。如果解析开始前创建的这个蛋白对象指针在函数返回前仍然为空，则其将在函数返回前被删除并析构
- Chain：解析开始时，创建新Model时，以及每次检测到链名发生变化时（从上一个"ATOM"行到当前行），都会创建一个新的链对象指针
- Residue：解析开始时，创建新Model时，创建新链时，以及每次检测到残基名、残基编号或残基插入字符三者之一发生变化时（从上一个"ATOM"行到当前行），都会创建一个新的残基对象指针
- Atom：每检测到一个新的"ATOM"行都会创建一个新的Atom对象指针
