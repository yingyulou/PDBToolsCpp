# PDBToolsCpp

PDB文件解析与坐标线性代数运算工具集。

PDB文件在PDBToolsCpp中将被解析为4个层级：Protein -> Chain -> Residue -> Atom

**本文档中所有的“角度”均指弧度制角度；所有的“旋转矩阵”均指右乘旋转矩阵。**

**任何结构对象均定义于堆内存上，使用时需注意内存管理。（详见下文析构函数部分内容）**

## 0. 编译及使用说明

* 导入"PDBTools"头文件即可使用：

``` Cpp
#include <PDBToolsCpp/PDBTools>
```

* 编译依赖：

1. boost
2. Eigen

* 编译器需支持GNU C++17或以上标准
* PDBToolsCpp的所有接口均位于namespace PDBTools下

## 1. PDB文件解析

### 1.1 load

``` Cpp
Protein *load(const string &pdbFilePath, bool parseHBool = false);
```

将PDB文件解析为Protein对象。

#### 参数：

* pdbFilePath：PDB文件路径
* parseHBool：是否开启氢原子解析

#### 返回值：

* Protein对象

#### 例：

``` Cpp
Protein *proPtr = load("xxxx.pdb");
```

### 1.2 loadModel

``` Cpp
vector<Protein *> loadModel(const string &pdbFilePath, bool parseHBool = false);
```

将含有"MODEL"关键词的PDB文件解析为Protein对象vector。

#### 参数：

* pdbFilePath：PDB文件路径
* parseHBool：是否开启氢原子解析

#### 返回值：

* Protein对象构成的vector

#### 例：

``` Cpp
vector<Protein *> proPtrList = loadModel("xxxx.pdb");
```

## 2. Protein

Protein类，用于表示一个蛋白。

### 2.1 Constructor

``` Cpp
explicit Protein(const string &name = "", int model = 0);
```

#### 参数：

* name：蛋白名
* model：Model编号

#### 例：

``` Cpp
auto proPtr = new Protein;
```

### 2.2 Getter / Setter

``` Cpp
string          &name ();
int              model();
vector<Chain *> &sub  ();

Protein *name (const string          &val);
Protein *model(int                    val);
Protein *sub  (const vector<Chain *> &val);
```

对应于Constructor各参数的Getter / Setter。

#### 例：

``` Cpp
auto proPtr = new Protein;

auto name  = proPtr->name();
auto model = proPtr->model();
auto sub   = proPtr->sub();

proPtr
    ->name ("")
    ->model(0)
    ->sub  ({});
```

### 2.3 copy

``` Cpp
Protein *copy();
```

得到this的深拷贝。

#### 参数：

* void

#### 返回值：

* this的深拷贝

#### 例：

``` Cpp
auto proPtr = new Protein;

auto copyProPtr = proPtr->copy();
```

### 2.4 getResidues

``` Cpp
vector<Residue *> getResidues();
```

得到this包含的所有残基。

#### 参数：

* void

#### 返回值：

* this包含的所有残基

#### 例：

``` Cpp
auto proPtr = new Protein;

auto resPtrList = proPtr->getResidues();
```

### 2.5 getAtoms

``` Cpp
vector<Atom *> getAtoms();
```

得到this包含的所有原子。

#### 参数：

* void

#### 返回值：

* this包含的所有原子

#### 例：

``` Cpp
auto proPtr = new Protein;

auto atomPtrList = proPtr->getAtoms();
```

### 2.6 subMap

``` Cpp
unordered_map<string, Chain *> subMap();
```

得到this包含的所有链名 -> 链对象哈希表。

#### 参数：

* void

#### 返回值：

* this包含的所有链名 -> 链对象哈希表

#### 例：

``` Cpp
auto proPtr = new Protein;

auto subMap = proPtr->subMap();
```

### 2.7 dump

``` Cpp
Protein *dump(const string &dumpFilePath, const string &fileMode = "w");
```

将this输出至PDB文件。

#### 参数：

* dumpFilePath：输出文件路径
* fileMode：文件打开模式

#### 返回值：

* this

#### 例：

``` Cpp
auto proPtr = new Protein;

proPtr->dump("xxx.pdb");
```

### 2.8 begin, end

``` Cpp
typename vector<Chain *>::iterator begin();
typename vector<Chain *>::iterator end();
```

委托至sub()的迭代器。

#### 例：

``` Cpp
auto proPtr = new Protein;

for (auto chainPtr: *proPtr);
```

### 2.9 filterAtoms

``` Cpp
vector<Atom *> filterAtoms(
    const unordered_set<string> &atomNameSet = {"CA"});
```

按原子名筛选this包含的所有原子。

#### 参数：

* atomNameSet：原子名集合

#### 返回值：

* 筛选出的所有原子

#### 例：

``` Cpp
auto proPtr = new Protein;

auto atomPtrList = proPtr->filterAtoms();
```

### 2.10 getAtomsCoord

``` Cpp
Matrix<double, Dynamic, 3> getAtomsCoord();
```

得到this包含的所有原子坐标。

#### 参数：

* void

#### 返回值：

* this包含的所有原子坐标

#### 例：

``` Cpp
auto proPtr = new Protein;

auto coordArray = proPtr->getAtomsCoord();
```

### 2.11 filterAtomsCoord

``` Cpp
Matrix<double, Dynamic, 3> filterAtomsCoord(
    const unordered_set<string> &atomNameSet = {"CA"});
```

按原子名筛选this包含的所有原子坐标。

#### 参数：

* atomNameSet：原子名集合

#### 返回值：

* 筛选出的所有原子坐标

#### 例：

``` Cpp
auto proPtr = new Protein;

auto coordArray = proPtr->filterAtomsCoord();
```

### 2.12 center

``` Cpp
RowVector3d center();
```

得到this包含的所有原子坐标的几何中心。

#### 参数：

* void

#### 返回值：

* this包含的所有原子坐标的几何中心

#### 例：

``` Cpp
auto proPtr = new Protein;

auto centerCoord = proPtr->center();
```

### 2.13 moveCenter

``` Cpp
Protein *moveCenter();
```

平移this的所有原子，使得其几何中心变为原点。

#### 参数：

* void

#### 返回值：

* this

#### 例：

``` Cpp
auto proPtr = new Protein;

proPtr->moveCenter();
```

### 2.14 seq

``` Cpp
string seq();
```

得到this的序列。

#### 参数：

* void

#### 返回值：

* this的序列

#### 例：

``` Cpp
auto proPtr = new Protein;

auto seqStr = proPtr->seq();
```

### 2.15 fastaStr

``` Cpp
string fastaStr(const string &titleStr = "");
```

得到this的Fasta格式字符串。

#### 参数：

* titleStr：Fasta标题。如果传入空字符串，则将自动为其分配一个标题

#### 返回值：

* this的Fasta格式字符串

#### 例：

``` Cpp
auto proPtr = new Protein;

auto fastaStr = proPtr->fastaStr();
```

### 2.16 dumpFasta

``` Cpp
Protein *dumpFasta(const string &dumpFilePath,
    const string &fileMode = "w", const string &titleStr = "");
```

将this输出至Fasta文件。

#### 参数：

* dumpFilePath：输出文件路径
* fileMode：文件打开模式
* titleStr：Fasta标题。如果传入空字符串，则将自动为其分配一个标题

#### 返回值：

* this

#### 例：

``` Cpp
auto proPtr = new Protein;

proPtr->dumpFasta("xxx.fasta");
```

### 2.17 renumResidues

``` Cpp
Protein *renumResidues(int startNum = 1);
```

对this包含的所有残基进行重编号。

#### 参数：

* startNum：起始编号

#### 返回值：

* this

#### 例：

``` Cpp
auto proPtr = new Protein;

proPtr->renumResidues();
```

### 2.18 renumAtoms

``` Cpp
Protein *renumAtoms(int startNum = 1);
```

对this包含的所有原子进行重编号。

#### 参数：

* startNum：起始编号

#### 返回值：

* this

#### 例：

``` Cpp
auto proPtr = new Protein;

proPtr->renumAtoms();
```

### 2.19 append

``` Cpp
Protein *append(Chain *subPtr, bool copyBool = true);
```

在sub()的末尾插入链对象。

#### 参数：

* subPtr：链对象
* copyBool：是否需要插入一个拷贝的链对象

#### 返回值：

* this

#### 例：

``` Cpp
auto proPtr   = new Protein;
auto chainPtr = new Chain;

proPtr->append(chainPtr);
```

### 2.20 insert

``` Cpp
Protein *insert(typename vector<Chain *>::iterator insertIter,
    Chain *subPtr, bool copyBool = true);
```

在sub()的任意位置插入链对象。

#### 参数：

* insertIter：插入位置迭代器
* subPtr：链对象
* copyBool：是否需要插入一个拷贝的链对象

#### 返回值：

* this

#### 例：

``` Cpp
auto proPtr   = new Protein;
auto chainPtr = new Chain;

proPtr->insert(proPtr->sub().begin(), chainPtr);
```

### 2.21 removeAlt

``` Cpp
Protein *removeAlt();
```

遍历this包含的所有原子对象，如果原子对象的alt()为""，则忽略；如果为"A"，则修改为""；否则，删除当前原子。

#### 参数：

* void

#### 返回值：

* this

#### 例：

``` Cpp
auto proPtr = new Protein;

proPtr->removeAlt();
```

### 2.22 dumpStr

``` Cpp
string dumpStr();
```

得到this的PDB格式的字符串。

#### 参数：

* void

#### 返回值：

* this的PDB格式的字符串

#### 例：

``` Cpp
auto proPtr = new Protein;

auto pdbStr = proPtr->dumpStr();
```

### 2.23 Destructor

``` Cpp
~Protein();
```

#### 例：

``` Cpp
auto proPtr = new Protein;

delete proPtr;
```

## 3. Chain

Chain类，用于表示一条链。

### 3.1 Constructor

``` Cpp
explicit Chain(const string &name = "", Protein *owner = nullptr);
```

#### 参数：

* name：链名
* owner：this所属的Protein

#### 例：

``` Cpp
auto chainPtr = new Chain;
```

### 3.2 Getter / Setter

``` Cpp
string            &name ();
Protein           *owner();
vector<Residue *> &sub  ();

Chain *name (const string            &val);
Chain *owner(Protein                 *val);
Chain *sub  (const vector<Residue *> &val);
```

对应于Constructor各参数的Getter / Setter。

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto name  = chainPtr->name();
auto owner = chainPtr->owner();
auto sub   = chainPtr->sub();

chainPtr
    ->name ("")
    ->owner(nullptr)
    ->sub  ({});
```

### 3.3 copy

``` Cpp
Chain *copy();
```

得到this的深拷贝。

#### 参数：

* void

#### 返回值：

* this的深拷贝

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto copyChainPtr = chainPtr->copy();
```

### 3.4 getResidues

``` Cpp
vector<Residue *> getResidues();
```

得到this包含的所有残基。

#### 参数：

* void

#### 返回值：

* this包含的所有残基

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto resPtrList = chainPtr->getResidues();
```

### 3.5 getAtoms

``` Cpp
vector<Atom *> getAtoms();
```

得到this包含的所有原子。

#### 参数：

* void

#### 返回值：

* this包含的所有原子

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto atomPtrList = chainPtr->getAtoms();
```

### 3.6 subMap

``` Cpp
unordered_map<string, Residue *> subMap();
```

得到this包含的所有残基完整编号 -> 残基对象哈希表。

#### 参数：

* void

#### 返回值：

* this包含的所有残基完整编号 -> 残基对象哈希表

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto subMap = chainPtr->subMap();
```

### 3.7 dump

``` Cpp
Chain *dump(const string &dumpFilePath, const string &fileMode = "w");
```

将this输出至PDB文件。

#### 参数：

* dumpFilePath：输出文件路径
* fileMode：文件打开模式

#### 返回值：

* this

#### 例：

``` Cpp
auto chainPtr = new Chain;

chainPtr->dump("xxx.pdb");
```

### 3.8 begin, end

``` Cpp
typename vector<Residue *>::iterator begin();
typename vector<Residue *>::iterator end();
```

委托至sub()的迭代器。

#### 例：

``` Cpp
auto chainPtr = new Chain;

for (auto resPtr: *chainPtr);
```

### 3.9 filterAtoms

``` Cpp
vector<Atom *> filterAtoms(
    const unordered_set<string> &atomNameSet = {"CA"});
```

按原子名筛选this包含的所有原子。

#### 参数：

* atomNameSet：原子名集合

#### 返回值：

* 筛选出的所有原子

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto atomPtrList = chainPtr->filterAtoms();
```

### 3.10 getAtomsCoord

``` Cpp
Matrix<double, Dynamic, 3> getAtomsCoord();
```

得到this包含的所有原子坐标。

#### 参数：

* void

#### 返回值：

* this包含的所有原子坐标

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto coordArray = chainPtr->getAtomsCoord();
```

### 3.11 filterAtomsCoord

``` Cpp
Matrix<double, Dynamic, 3> filterAtomsCoord(
    const unordered_set<string> &atomNameSet = {"CA"});
```

按原子名筛选this包含的所有原子坐标。

#### 参数：

* atomNameSet：原子名集合

#### 返回值：

* 筛选出的所有原子坐标

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto coordArray = chainPtr->filterAtomsCoord();
```

### 3.12 center

``` Cpp
RowVector3d center();
```

得到this包含的所有原子坐标的几何中心。

#### 参数：

* void

#### 返回值：

* this包含的所有原子坐标的几何中心

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto centerCoord = chainPtr->center();
```

### 3.13 moveCenter

``` Cpp
Chain *moveCenter();
```

平移this的所有原子，使得其几何中心变为原点。

#### 参数：

* void

#### 返回值：

* this

#### 例：

``` Cpp
auto chainPtr = new Chain;

chainPtr->moveCenter();
```

### 3.14 seq

``` Cpp
string seq();
```

得到this的序列。

#### 参数：

* void

#### 返回值：

* this的序列

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto seqStr = chainPtr->seq();
```

### 3.15 fastaStr

``` Cpp
string fastaStr(const string &titleStr = "");
```

得到this的Fasta格式字符串。

#### 参数：

* titleStr：Fasta标题。如果传入空字符串，则将自动为其分配一个标题

#### 返回值：

* this的Fasta格式字符串

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto fastaStr = chainPtr->fastaStr();
```

### 3.16 dumpFasta

``` Cpp
Chain *dumpFasta(const string &dumpFilePath,
    const string &fileMode = "w", const string &titleStr = "");
```

将this输出至Fasta文件。

#### 参数：

* dumpFilePath：输出文件路径
* fileMode：文件打开模式
* titleStr：Fasta标题。如果传入空字符串，则将自动为其分配一个标题

#### 返回值：

* this

#### 例：

``` Cpp
auto chainPtr = new Chain;

chainPtr->dumpFasta("xxx.fasta");
```

### 3.17 renumResidues

``` Cpp
Chain *renumResidues(int startNum = 1);
```

对this包含的所有残基进行重编号。

#### 参数：

* startNum：起始编号

#### 返回值：

* this

#### 例：

``` Cpp
auto chainPtr = new Chain;

chainPtr->renumResidues();
```

### 3.18 renumAtoms

``` Cpp
Chain *renumAtoms(int startNum = 1);
```

对this包含的所有原子进行重编号。

#### 参数：

* startNum：起始编号

#### 返回值：

* this

#### 例：

``` Cpp
auto chainPtr = new Chain;

chainPtr->renumAtoms();
```

### 3.19 append

``` Cpp
Chain *append(Residue *subPtr, bool copyBool = true);
```

在sub()的末尾插入残基对象。

#### 参数：

* subPtr：残基对象
* copyBool：是否需要插入一个拷贝的残基对象

#### 返回值：

* this

#### 例：

``` Cpp
auto chainPtr = new Chain;
auto resPtr   = new Residue;

chainPtr->append(resPtr);
```

### 3.20 insert

``` Cpp
Chain *insert(typename vector<Residue *>::iterator insertIter,
    Residue *subPtr, bool copyBool = true);
```

在sub()的任意位置插入链对象。

#### 参数：

* insertIter：插入位置迭代器
* subPtr：残基对象
* copyBool：是否需要插入一个拷贝的残基对象

#### 返回值：

* this

#### 例：

``` Cpp
auto chainPtr = new Chain;
auto resPtr   = new Residue;

chainPtr->insert(chainPtr->sub().begin(), resPtr);
```

### 3.21 removeAlt

``` Cpp
Chain *removeAlt();
```

遍历this包含的所有原子对象，如果原子对象的alt()为""，则忽略；如果为"A"，则修改为""；否则，删除当前原子。

#### 参数：

* void

#### 返回值：

* this

#### 例：

``` Cpp
auto chainPtr = new Chain;

chainPtr->removeAlt();
```

### 3.22 dumpStr

``` Cpp
string dumpStr();
```

得到this的PDB格式的字符串。

#### 参数：

* void

#### 返回值：

* this的PDB格式的字符串

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto pdbStr = chainPtr->dumpStr();
```

### 3.23 iter

``` Cpp
typename vector<Chain *>::iterator iter();
```

得到this在this->owner()->sub()中的迭代器。

#### 参数：

* void

#### 返回值：

* this在this->owner()->sub()中的迭代器

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto chainIter = chainPtr->iter();
```

### 3.24 prev

``` Cpp
Chain *prev(int shiftLen = 1);
```

得到this在this->owner()->sub()中的向前第N个链对象。

#### 参数：

* shiftLen：向前偏移量

#### 返回值：

* this在this->owner()->sub()中的向前第N个链对象

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto prevChainPtr = chainPtr->prev();
```

### 3.25 next

``` Cpp
Chain *next(int shiftLen = 1);
```

得到this在this->owner()->sub()中的向后第N个链对象。

#### 参数：

* shiftLen：向后偏移量

#### 返回值：

* this在this->owner()->sub()中的向后第N个链对象

#### 例：

``` Cpp
auto chainPtr = new Chain;

auto nextChainPtr = chainPtr->next();
```

### 3.26 remove

``` Cpp
typename vector<Chain *>::iterator remove(bool deteleBool = true);
```

从this->owner()->sub()中删除this。

#### 参数：

* deteleBool：是否需要析构this

#### 返回值：

* this的后继元素迭代器

#### 例：

``` Cpp
auto chainPtr = new Chain;

chainPtr->remove();
```

### 3.27 Destructor

``` Cpp
~Chain();
```

#### 例：

``` Cpp
auto chainPtr = new Chain;

delete chainPtr;
```

## 4. Residue

Residue类，用于表示一个残基。

### 4.1 Constructor

``` Cpp
explicit Residue(const string &name = "", int num = 0,
    const string &ins = "", Chain *owner = nullptr);
```

#### 参数：

* name：残基名
* num：残基编号
* ins：残基插入字符
* owner：this所属的Chain

#### 例：

``` Cpp
auto resPtr = new Residue;
```

### 4.2 Getter / Setter

``` Cpp
string         &name ();
int             num  ();
string         &ins  ();
Chain          *owner();
vector<Atom *> &sub  ();

Residue *name (const string         &val);
Residue *num  (int                   val);
Residue *ins  (const string         &val);
Residue *owner(Chain                *val);
Residue *sub  (const vector<Atom *> &val);
```

对应于Constructor各参数的Getter / Setter。

#### 例：

``` Cpp
auto resPtr = new Residue;

auto name  = resPtr->name();
auto num   = resPtr->num();
auto ins   = resPtr->ins();
auto owner = resPtr->owner();
auto sub   = resPtr->sub();

resPtr
    ->name ("")
    ->num  (0)
    ->ins  ("")
    ->owner(nullptr)
    ->sub  ({});
```

### 4.3 compNum

``` Cpp
string compNum();

Residue *compNum(int num, const string &ins = "");
Residue *compNum(const pair<int, string> &compNumPair);
```

直接获取或设定残基完整编号。

#### 参数：

* num：残基编号
* ins：残基插入字符
* compNumPair：残基完整编号

#### 返回值：

* 残基完整编号 / this

#### 例：

``` Cpp
auto resPtr = new Residue;

auto compNum = resPtr->compNum();

resPtr
    ->compNum(0, "")
    ->compNum({0, ""});
```

### 4.4 copy

``` Cpp
Residue *copy();
```

得到this的深拷贝。

#### 参数：

* void

#### 返回值：

* this的深拷贝

#### 例：

``` Cpp
auto resPtr = new Residue;

auto copyResPtr = resPtr->copy();
```

### 4.5 getResidues

``` Cpp
vector<Residue *> getResidues();
```

得到this包含的所有残基。

#### 参数：

* void

#### 返回值：

* this包含的所有残基

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtrList = resPtr->getResidues();
```

### 4.6 getAtoms

``` Cpp
vector<Atom *> getAtoms();
```

得到this包含的所有原子。

#### 参数：

* void

#### 返回值：

* this包含的所有原子

#### 例：

``` Cpp
auto resPtr = new Residue;

auto atomPtrList = resPtr->getAtoms();
```

### 4.7 subMap

``` Cpp
unordered_map<string, Atom *> subMap();
```

得到this包含的所有原子名 -> 原子对象哈希表。

#### 参数：

* void

#### 返回值：

* this包含的所有原子名 -> 原子对象哈希表

#### 例：

``` Cpp
auto resPtr = new Residue;

auto subMap = resPtr->subMap();
```

### 4.8 coordMap

``` Cpp
unordered_map<string, RowVector3d> coordMap();
```

得到this包含的所有原子名 -> 原子坐标哈希表。

#### 参数：

* void

#### 返回值：

* this包含的所有原子名 -> 原子坐标哈希表

#### 例：

``` Cpp
auto resPtr = new Residue;

auto coordMap = resPtr->coordMap();
```

### 4.9 calcBBDihedralAngle

``` Cpp
double calcBBDihedralAngle(DIH dihedralEnum);
```

计算主链二面角。

#### 参数：

* dihedralEnum：主链二面角种类。DIH::PHI或DIH::L表示Phi；DIH::PSI或DIH::R表示Psi

#### 返回值：

* 有符号二面角值（-pi ~ pi）

#### 例：

``` Cpp
auto resPtr = new Residue;

auto dihedralAngle = resPtr->calcBBDihedralAngle(DIH::L);
```

### 4.10 calcBBRotationMatrixByDeltaAngle

``` Cpp
pair<RowVector3d, Matrix3d> calcBBRotationMatrixByDeltaAngle(
    DIH dihedralEnum, SIDE sideEnum, double deltaAngle);
```

以旋转角度作为参数，计算主链旋转矩阵。

#### 参数：

* dihedralEnum：主链二面角种类。DIH::PHI或DIH::L表示Phi；DIH::PSI或DIH::R表示Psi
* sideEnum：转动侧。SIDE::N或SIDE::L表示转动N端；SIDE::C或SIDE::R表示转动C端
* deltaAngle：旋转角度

#### 返回值：

* 旋转前/后平移向量
* 旋转矩阵

#### 例：

``` Cpp
auto resPtr = new Residue;

auto [moveCoord, rotationMatrix] = resPtr->calcBBRotationMatrixByDeltaAngle(DIH::PHI, SIDE::N, 1.);
```

### 4.11 calcBBRotationMatrixByTargetAngle

``` Cpp
pair<RowVector3d, Matrix3d> calcBBRotationMatrixByTargetAngle(
    DIH dihedralEnum, SIDE sideEnum, double targetAngle);
```

以目标角度作为参数，计算主链旋转矩阵。

#### 参数：

* dihedralEnum：主链二面角种类。DIH::PHI或DIH::L表示Phi；DIH::PSI或DIH::R表示Psi
* sideEnum：转动侧。SIDE::N或SIDE::L表示转动N端；SIDE::C或SIDE::R表示转动C端
* targetAngle：目标角度

#### 返回值：

* 旋转前/后平移向量
* 旋转矩阵

#### 例：

``` Cpp
auto resPtr = new Residue;

auto [moveCoord, rotationMatrix] = resPtr->calcBBRotationMatrixByTargetAngle(DIH::PHI, SIDE::N, 0.);
```

### 4.12 getBBRotationAtomPtr

``` Cpp
vector<Atom *> getBBRotationAtomPtr(DIH dihedralEnum, SIDE sideEnum);
```

获取以给定参数进行旋转时，所有需要旋转的原子对象列表。

#### 参数：

* dihedralEnum：主链二面角种类。DIH::PHI或DIH::L表示Phi；DIH::PSI或DIH::R表示Psi
* sideEnum：转动侧。SIDE::N或SIDE::L表示转动N端；SIDE::C或SIDE::R表示转动C端

#### 返回值：

* 所有需要旋转的原子对象列表

#### 例：

``` Cpp
auto resPtr = new Residue;

auto rotationAtomPtrList = resPtr->getBBRotationAtomPtr(DIH::PHI, SIDE::N);
```

### 4.13 rotateBBDihedralAngleByDeltaAngle

``` Cpp
Residue *rotateBBDihedralAngleByDeltaAngle(DIH dihedralEnum,
    SIDE sideEnum, double deltaAngle);
```

以旋转角度作为参数直接旋转主链。

#### 参数：

* dihedralEnum：主链二面角种类。DIH::PHI或DIH::L表示Phi；DIH::PSI或DIH::R表示Psi
* sideEnum：转动侧。SIDE::N或SIDE::L表示转动N端；SIDE::C或SIDE::R表示转动C端
* deltaAngle：旋转角度

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtr->rotateBBDihedralAngleByDeltaAngle(DIH::PHI, SIDE::N, 1.);
```

### 4.14 rotateBBDihedralAngleByTargetAngle

``` Cpp
Residue *rotateBBDihedralAngleByTargetAngle(DIH dihedralEnum,
    SIDE sideEnum, double targetAngle);
```

以目标角度作为参数直接旋转主链。

#### 参数：

* dihedralEnum：主链二面角种类。DIH::PHI或DIH::L表示Phi；DIH::PSI或DIH::R表示Psi
* sideEnum：转动侧。SIDE::N或SIDE::L表示转动N端；SIDE::C或SIDE::R表示转动C端
* targetAngle：目标角度

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtr->rotateBBDihedralAngleByTargetAngle(DIH::PHI, SIDE::N, 0.);
```

### 4.15 calcSCDihedralAngle

``` Cpp
double calcSCDihedralAngle(int dihedralIdx);
```

计算侧链二面角。

#### 参数：

* dihedralIdx：侧链二面角索引值。索引值从0开始编号，最大允许索引值根据残基种类而不同。索引值表示某个残基从主链到侧链方向上的第N个侧链二面角

#### 返回值：

* 有符号二面角值（-pi ~ pi）

#### 例：

``` Cpp
auto resPtr = new Residue;

double dihedralAngle = resPtr->calcSCDihedralAngle(0);
```

### 4.16 calcSCRotationMatrixByDeltaAngle

``` Cpp
pair<RowVector3d, Matrix3d> calcSCRotationMatrixByDeltaAngle(
    int dihedralIdx, double deltaAngle);
```

以旋转角度/目标角度作为参数，计算侧链旋转矩阵。

#### 参数：

* dihedralIdx：侧链二面角索引值。索引值从0开始编号，最大允许索引值根据残基种类而不同。索引值表示某个残基从主链到侧链方向上的第N个侧链二面角
* deltaAngle：旋转角度

#### 返回值：

* 旋转前/后平移向量
* 旋转矩阵

#### 例：

``` Cpp
auto resPtr = new Residue;

auto [moveCoord, rotationMatrix] = resPtr->calcSCRotationMatrixByDeltaAngle(0, 1.);
```

### 4.17 calcSCRotationMatrixByTargetAngle

``` Cpp
pair<RowVector3d, Matrix3d> calcSCRotationMatrixByTargetAngle(
    int dihedralIdx, double targetAngle);
```

以目标角度作为参数，计算侧链旋转矩阵。

#### 参数：

* dihedralIdx：侧链二面角索引值。索引值从0开始编号，最大允许索引值根据残基种类而不同。索引值表示某个残基从主链到侧链方向上的第N个侧链二面角
* targetAngle：目标角度

#### 返回值：

* 旋转前/后平移向量
* 旋转矩阵

#### 例：

``` Cpp
auto resPtr = new Residue;

auto [moveCoord, rotationMatrix] = resPtr->calcSCRotationMatrixByTargetAngle(0, 0.);
```

### 4.18 getSCRotationAtomPtr

``` Cpp
vector<Atom *> getSCRotationAtomPtr(int dihedralIdx);
```

获取以给定侧链二面角进行旋转时，所有需要旋转的原子对象列表。

#### 参数：

* dihedralIdx：侧链二面角索引值。索引值从0开始编号，最大允许索引值根据残基种类而不同。索引值表示某个残基从主链到侧链方向上的第N个侧链二面角

#### 返回值：

* 所有需要旋转的原子对象列表

#### 例：

``` Cpp
auto resPtr = new Residue;

auto rotationAtomPtrList = resPtr->getSCRotationAtomPtr(0);
```

### 4.19 rotateSCDihedralAngleByDeltaAngle

``` Cpp
Residue *rotateSCDihedralAngleByDeltaAngle(int dihedralIdx, double deltaAngle);
```

以旋转角度作为参数直接旋转侧链。

#### 参数：

* dihedralIdx：侧链二面角索引值。索引值从0开始编号，最大允许索引值根据残基种类而不同。索引值表示某个残基从主链到侧链方向上的第N个侧链二面角
* deltaAngle：旋转角度

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtr->rotateSCDihedralAngleByDeltaAngle(0, 1.);
```

### 4.20 rotateSCDihedralAngleByTargetAngle

``` Cpp
Residue *rotateSCDihedralAngleByTargetAngle(int dihedralIdx, double targetAngle);
```

以目标角度作为参数直接旋转侧链。

#### 参数：

* dihedralIdx：侧链二面角索引值。索引值从0开始编号，最大允许索引值根据残基种类而不同。索引值表示某个残基从主链到侧链方向上的第N个侧链二面角
* targetAngle：目标角度

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtr->rotateSCDihedralAngleByTargetAngle(0, 0.);
```

### 4.21 dump

``` Cpp
Residue *dump(const string &dumpFilePath, const string &fileMode = "w");
```

将this输出至PDB文件。

#### 参数：

* dumpFilePath：输出文件路径
* fileMode：文件打开模式

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtr->dump("xxx.pdb");
```

### 4.22 begin, end

``` Cpp
typename vector<Atom *>::iterator begin();
typename vector<Atom *>::iterator end();
```

委托至sub()的迭代器。

#### 例：

``` Cpp
auto resPtr = new Residue;

for (auto atomPtr: *resPtr);
```

### 4.23 filterAtoms

``` Cpp
vector<Atom *> filterAtoms(
    const unordered_set<string> &atomNameSet = {"CA"});
```

按原子名筛选this包含的所有原子。

#### 参数：

* atomNameSet：原子名集合

#### 返回值：

* 筛选出的所有原子

#### 例：

``` Cpp
auto resPtr = new Residue;

auto atomPtrList = resPtr->filterAtoms();
```

### 4.24 getAtomsCoord

``` Cpp
Matrix<double, Dynamic, 3> getAtomsCoord();
```

得到this包含的所有原子坐标。

#### 参数：

* void

#### 返回值：

* this包含的所有原子坐标

#### 例：

``` Cpp
auto resPtr = new Residue;

auto coordArray = resPtr->getAtomsCoord();
```

### 4.25 filterAtomsCoord

``` Cpp
Matrix<double, Dynamic, 3> filterAtomsCoord(
    const unordered_set<string> &atomNameSet = {"CA"});
```

按原子名筛选this包含的所有原子坐标。

#### 参数：

* atomNameSet：原子名集合

#### 返回值：

* 筛选出的所有原子坐标

#### 例：

``` Cpp
auto resPtr = new Residue;

auto coordArray = resPtr->filterAtomsCoord();
```

### 4.26 center

``` Cpp
RowVector3d center();
```

得到this包含的所有原子坐标的几何中心。

#### 参数：

* void

#### 返回值：

* this包含的所有原子坐标的几何中心

#### 例：

``` Cpp
auto resPtr = new Residue;

auto centerCoord = resPtr->center();
```

### 4.27 moveCenter

``` Cpp
Residue *moveCenter();
```

平移this的所有原子，使得其几何中心变为原点。

#### 参数：

* void

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtr->moveCenter();
```

### 4.28 seq

``` Cpp
string seq();
```

得到this的序列。

#### 参数：

* void

#### 返回值：

* this的序列

#### 例：

``` Cpp
auto resPtr = new Residue;

auto seqStr = resPtr->seq();
```

### 4.29 fastaStr

``` Cpp
string fastaStr(const string &titleStr = "");
```

得到this的Fasta格式字符串。

#### 参数：

* titleStr：Fasta标题。如果传入空字符串，则将自动为其分配一个标题

#### 返回值：

* this的Fasta格式字符串

#### 例：

``` Cpp
auto resPtr = new Residue;

auto fastaStr = resPtr->fastaStr();
```

### 4.30 dumpFasta

``` Cpp
Residue *dumpFasta(const string &dumpFilePath,
    const string &fileMode = "w", const string &titleStr = "");
```

将this输出至Fasta文件。

#### 参数：

* dumpFilePath：输出文件路径
* fileMode：文件打开模式
* titleStr：Fasta标题。如果传入空字符串，则将自动为其分配一个标题

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtr->dumpFasta("xxx.fasta");
```

### 4.31 renumResidues

``` Cpp
Residue *renumResidues(int startNum = 1);
```

对this包含的所有残基进行重编号。

#### 参数：

* startNum：起始编号

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtr->renumResidues();
```

### 4.32 renumAtoms

``` Cpp
Residue *renumAtoms(int startNum = 1);
```

对this包含的所有原子进行重编号。

#### 参数：

* startNum：起始编号

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtr->renumAtoms();
```

### 4.33 append

``` Cpp
Residue *append(Atom *subPtr, bool copyBool = true);
```

在sub()的末尾插入原子对象。

#### 参数：

* subPtr：原子对象
* copyBool：是否需要插入一个拷贝的原子对象

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr  = new Residue;
auto atomPtr = new Atom;

resPtr->append(atomPtr);
```

### 4.34 insert

``` Cpp
Residue *insert(typename vector<Atom *>::iterator insertIter,
    Atom *subPtr, bool copyBool = true);
```

在sub()的任意位置插入原子对象。

#### 参数：

* insertIter：插入位置迭代器
* subPtr：原子对象
* copyBool：是否需要插入一个拷贝的原子对象

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr  = new Residue;
auto atomPtr = new Atom;

resPtr->insert(resPtr->sub().begin(), atomPtr);
```

### 4.35 removeAlt

``` Cpp
Residue *removeAlt();
```

遍历this包含的所有原子对象，如果原子对象的alt()为""，则忽略；如果为"A"，则修改为""；否则，删除当前原子。

#### 参数：

* void

#### 返回值：

* this

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtr->removeAlt();
```

### 4.36 dumpStr

``` Cpp
string dumpStr();
```

得到this的PDB格式的字符串。

#### 参数：

* void

#### 返回值：

* this的PDB格式的字符串

#### 例：

``` Cpp
auto resPtr = new Residue;

auto pdbStr = resPtr->dumpStr();
```

### 4.37 iter

``` Cpp
typename vector<Residue *>::iterator iter();
```

得到this在this->owner()->sub()中的迭代器。

#### 参数：

* void

#### 返回值：

* this在this->owner()->sub()中的迭代器

#### 例：

``` Cpp
auto resPtr = new Residue;

auto resIter = resPtr->iter();
```

### 4.38 prev

``` Cpp
Residue *prev(int shiftLen = 1);
```

得到this在this->owner()->sub()中的向前第N个残基对象。

#### 参数：

* shiftLen：向前偏移量

#### 返回值：

* this在this->owner()->sub()中的向前第N个残基对象

#### 例：

``` Cpp
auto resPtr = new Residue;

auto prevResPtr = resPtr->prev();
```

### 4.39 next

``` Cpp
Residue *next(int shiftLen = 1);
```

得到this在this->owner()->sub()中的向后第N个残基对象。

#### 参数：

* shiftLen：向后偏移量

#### 返回值：

* this在this->owner()->sub()中的向后第N个残基对象

#### 例：

``` Cpp
auto resPtr = new Residue;

auto nextResPtr = resPtr->next();
```

### 4.40 remove

``` Cpp
typename vector<Residue *>::iterator remove(bool deteleBool = true);
```

从this->owner()->sub()中删除this。

#### 参数：

* deteleBool：是否需要析构this

#### 返回值：

* this的后继元素迭代器

#### 例：

``` Cpp
auto resPtr = new Residue;

resPtr->remove();
```

### 4.41 Destructor

``` Cpp
~Residue();
```

#### 例：

``` Cpp
auto resPtr = new Residue;

delete resPtr;
```

## 5. Atom

Atom类，用于表示一个原子。

### 5.1 Constructor

``` Cpp
explicit Atom(const string &name = "", int num = 0,
    const RowVector3d &coord = RowVector3d::Zero(),
    const string &alt = "", const string &occ = "",
    const string &tempF = "", const string &ele = "",
    const string &chg = "", Residue *owner = nullptr);
```

#### 参数：

* name：原子名
* num：原子编号
* coord：原子坐标
* alt：备用位置指示符
* occ：占有
* tempF：温度因子
* ele：元素符号
* chg：电荷
* owner：this所属的Residue

#### 例：

``` Cpp
auto atomPtr = new Atom;
```

### 5.2 Getter / Setter

``` Cpp
string      &name ();
int          num  ();
RowVector3d &coord();
string      &alt  ();
string      &occ  ();
string      &tempF();
string      &ele  ();
string      &chg  ();
Residue     *owner();

Atom *name (const string      &val);
Atom *num  (int                val);
Atom *coord(const RowVector3d &val);
Atom *alt  (const string      &val);
Atom *occ  (const string      &val);
Atom *tempF(const string      &val);
Atom *ele  (const string      &val);
Atom *chg  (const string      &val);
Atom *owner(Residue           *val);
```

对应于Constructor各参数的Getter / Setter。

#### 例：

``` Cpp
auto atomPtr = new Atom;

auto name  = atomPtr->name ();
auto num   = atomPtr->num  ();
auto coord = atomPtr->coord();
auto alt   = atomPtr->alt  ();
auto occ   = atomPtr->occ  ();
auto tempF = atomPtr->tempF();
auto ele   = atomPtr->ele  ();
auto chg   = atomPtr->chg  ();
auto owner = atomPtr->owner();

atomPtr
    ->name ("")
    ->num  (0)
    ->coord(RowVector3d::Zero())
    ->alt  ("")
    ->occ  ("")
    ->tempF("")
    ->ele  ("")
    ->chg  ("")
    ->owner(nullptr);
```

### 5.3 copy

``` Cpp
Atom *copy();
```

得到this的深拷贝。

#### 参数：

* void

#### 返回值：

* this的深拷贝

#### 例：

``` Cpp
auto atomPtr = new Atom;

auto copyAtomPtr = atomPtr->copy();
```

### 5.4 operator-

``` Cpp
double operator-(const Atom &rhs) const;
```

计算两原子间的欧几里得距离。

#### 参数：

* rhs：另一个原子对象

#### 返回值：

* 两原子间的欧几里得距离

#### 例：

``` Cpp
auto atomPtr = new Atom;

cout << *atomPtr - *atomPtr;
```

### 5.5 dump

``` Cpp
Atom *dump(const string &dumpFilePath, const string &fileMode = "w");
```

将this输出至PDB文件。

#### 参数：

* dumpFilePath：输出文件路径
* fileMode：文件打开模式

#### 返回值：

* this

#### 例：

``` Cpp
auto atomPtr = new Atom;

atomPtr->dump("xxx.pdb");
```

### 5.6 dumpStr

``` Cpp
string dumpStr();
```

得到this的PDB格式的字符串。

#### 参数：

* void

#### 返回值：

* this的PDB格式的字符串

#### 例：

``` Cpp
auto atomPtr = new Atom;

auto pdbStr = atomPtr->dumpStr();
```

### 5.7 iter

``` Cpp
typename vector<Atom *>::iterator iter();
```

得到this在this->owner()->sub()中的迭代器。

#### 参数：

* void

#### 返回值：

* this在this->owner()->sub()中的迭代器

#### 例：

``` Cpp
auto atomPtr = new Atom;

auto atomIter = atomPtr->iter();
```

### 5.8 prev

``` Cpp
Atom *prev(int shiftLen = 1);
```

得到this在this->owner()->sub()中的向前第N个残基对象。

#### 参数：

* shiftLen：向前偏移量

#### 返回值：

* this在this->owner()->sub()中的向前第N个残基对象

#### 例：

``` Cpp
auto atomPtr = new Atom;

auto prevAtomPtr = atomPtr->prev();
```

### 5.9 next

``` Cpp
Atom *next(int shiftLen = 1);
```

得到this在this->owner()->sub()中的向后第N个残基对象。

#### 参数：

* shiftLen：向后偏移量

#### 返回值：

* this在this->owner()->sub()中的向后第N个残基对象

#### 例：

``` Cpp
auto atomPtr = new Atom;

auto nextAtomPtr = atomPtr->next();
```

### 5.10 remove

``` Cpp
typename vector<Atom *>::iterator remove(bool deteleBool = true);
```

从this->owner()->sub()中删除this。

#### 参数：

* deteleBool：是否需要析构this

#### 返回值：

* this的后继元素迭代器

#### 例：

``` Cpp
auto atomPtr = new Atom;

atomPtr->remove();
```

## 6. 数学函数

### 6.1 degrees, radians

``` Cpp
double degrees(double radiansAngle);
double radians(double degreesAngle);
```

角度，弧度互相转换。

#### 参数：

* radiansAngle，degreesAngle：弧度值/角度值

#### 返回值：

* 转换后的角度值/弧度值

#### 例：

``` Cpp
cout << degrees(1.) << endl;
cout << radians(1.) << endl;
```

### 6.2 calcVectorAngle

``` Cpp
double calcVectorAngle(const RowVector3d &coordA, const RowVector3d &coordB);
```

计算两向量夹角。

#### 参数：

* coordA，coordB：两个三维向量

#### 返回值：

* 两向量夹角值（0 ~ pi）

#### 例：

``` Cpp
auto vectorAngle = calcVectorAngle(RowVector3d(1., 2., 3.), RowVector3d(4., 5., 6.));
```

### 6.3 calcRotationMatrix

``` Cpp
Matrix3d calcRotationMatrix(const RowVector3d &rotationAxis, double rotationAngle);
```

计算轴角旋转矩阵。

#### 参数：

* rotationAxis：旋转轴向量，无需缩放至单位长度
* rotationAngle：旋转角

#### 返回值：

* 旋转矩阵

#### 例：

``` Cpp
auto rotationMatrix = calcRotationMatrix(RowVector3d(1., 2., 3.), 1.);
```

### 6.4 calcRotationMatrixByTwoVector

``` Cpp
Matrix3d calcRotationMatrixByTwoVector(const RowVector3d &refCoord,
    const RowVector3d &tarCoord);
```

计算从向量tarCoord旋转至向量refCoord所需要的旋转矩阵。

#### 参数：

* refCoord, tarCoord：两个三维向量

#### 返回值：

* 旋转矩阵

#### 例：

``` Cpp
auto rotationMatrix = calcRotationMatrixByTwoVector(
    RowVector3d(1., 2., 3.), RowVector3d(4., 5., 6.));
```

### 6.5 calcDihedralAngle

``` Cpp
double calcDihedralAngle(
    const RowVector3d &coordA, const RowVector3d &coordB,
    const RowVector3d &coordC, const RowVector3d &coordD);
```

计算二面角。

#### 参数：

* coordA，coordB，coordC，coordD：四个三维向量

#### 返回值：

* 有符号二面角值（-pi ~ pi）

#### 例：

``` Cpp
auto dihedralAngle = calcDihedralAngle(
    RowVector3d(1., 2., 3.), RowVector3d(4., 5., 6.),
    RowVector3d(7., 8., 9.), RowVector3d(10., 11., 12.));
```

### 6.6 calcRMSD

``` Cpp
double calcRMSD(const Matrix<double, Dynamic, 3> &coordArrayA,
    const Matrix<double, Dynamic, 3> &coordArrayB);
```

对两组等长的三维坐标计算RMSD。

#### 参数：

* coordArrayA，coordArrayB：两组等长的矩阵（N * 3）

#### 返回值：

* RMSD值

#### 例：

``` Cpp
Matrix<double, 2, 3> coordArrayA, coordArrayB;

coordArrayA << 1., 2., 3., 4., 5., 6.;
coordArrayB << 7., 8., 9., 10., 11., 12.;

auto rmsdValue = calcRMSD(coordArrayA, coordArrayB);
```

### 6.7 calcSuperimposeRotationMatrix

``` Cpp
tuple<RowVector3d, Matrix3d, RowVector3d> calcSuperimposeRotationMatrix(
    const Matrix<double, Dynamic, 3> &tarCoordArray,
    const Matrix<double, Dynamic, 3> &srcCoordArray)
```

计算从srcCoordArray到tarCoordArray的叠合旋转矩阵。

此函数将使得((srcCoordArray.rowwise() - 平移向量A) * 旋转矩阵).rowwise() + 平移向量B与tarCoordArray形成叠合（RMSD最小）。

#### 参数：

* tarCoordArray, srcCoordArray：两组等长的矩阵（N * 3）

#### 返回值：

* 平移向量A
* 旋转矩阵
* 平移向量B

#### 例：

``` Cpp
Matrix<double, 2, 3> srcCoordArray, tarCoordArray;

srcCoordArray << 1., 2., 3., 4., 5., 6.;
tarCoordArray << 7., 8., 9., 10., 11., 12.;

auto [srcCenterCoord, rotationMatrix, tarCenterCoord] =
    calcSuperimposeRotationMatrix(tarCoordArray, srcCoordArray);

cout << ((srcCoordArray.rowwise() - srcCenterCoord) * rotationMatrix).rowwise() +
    tarCenterCoord << endl << tarCoordArray << endl;
```

### 6.8 calcRMSDAfterSuperimpose

``` Cpp
double calcRMSDAfterSuperimpose(
    const Matrix<double, Dynamic, 3> &tarCoordArray,
    const Matrix<double, Dynamic, 3> &srcCoordArray)
```

叠合并计算RMSD。

此函数会将srcCoordArray通过calcSuperimposeRotationMatrix函数向tarCoordArray进行叠合，然后计算两组坐标之间的RMSD。

#### 参数：

* tarCoordArray, srcCoordArray：两组等长的矩阵（N * 3）

#### 返回值：

* RMSD值

#### 例：

``` Cpp
Matrix<double, 2, 3> srcCoordArray, tarCoordArray;

srcCoordArray << 1., 2., 3., 4., 5., 6.;
tarCoordArray << 7., 8., 9., 10., 11., 12.;

auto rmsdValue = calcRMSDAfterSuperimpose(tarCoordArray, srcCoordArray);
```

## 7. 其他函数

### 7.1 operator<<

``` Cpp
ostream &operator<<(ostream &os, const Protein &proObj);
ostream &operator<<(ostream &os, const Chain   &chainObj);
ostream &operator<<(ostream &os, const Residue &resObj);
ostream &operator<<(ostream &os, const Atom    &atomObj);
```

任何结构对象均可用通过operator\<\<输出此对象的信息摘要。

#### 参数：

* os：输出流对象
* proObj, chainObj, resObj, atomObj：任何层级对象

#### 返回值：

* 输出流对象

#### 例：

``` Cpp
auto proPtr = new Protein;

cout << *proPtr << endl;
```

### 7.2 isH

``` Cpp
bool isH(const string &atomName);
```

判断一个原子名是否为氢原子。

#### 参数：

* atomName：原子名

#### 返回值：

* 原子名是否为氢原子

#### 例：

``` Cpp
bool isHBool = isH("1H");
```

### 7.3 splitCompNum

``` Cpp
pair<int, string> splitCompNum(const string &compNumStr);
```

将完整残基编号分割为残基编号和残基插入编号。

#### 参数：

* compNumStr：完整残基编号

#### 返回值：

* 残基编号
* 残基插入编号

#### 例：

``` Cpp
auto [resNum, resIns] = splitCompNum("1A");
```

### 7.4 dumpStr

``` Cpp
template <typename T>
string dumpStr(const T &structPtrList);
```

得到任意结构对象列表的PDB格式的字符串。

#### 参数：

* structPtrList：任意结构对象列表

#### 返回值：

* 任意结构对象列表的PDB格式的字符串

#### 例：

``` Cpp
auto proPtr = new Protein;

auto pdbStr = dumpStr(proPtr->sub());
```

### 7.5 dump

``` Cpp
template <typename T>
void dump(const T &structPtrList, const string &dumpFilePath,
    const string &fileMode = "w")
```

将任意结构对象列表输出至PDB文件。

#### 参数：

* structPtrList：任意结构对象列表
* dumpFilePath：输出文件路径
* fileMode：文件打开模式

#### 返回值：

* void

#### 例：

``` Cpp
auto proPtr = new Protein;

dump(proPtr->sub(), "xxx.pdb");
```

### 7.6 fastaStr

``` Cpp
template <typename T>
string fastaStr(const T &structPtrList)
```

得到任意结构对象列表的Fasta格式字符串。

#### 参数：

* structPtrList：任意结构对象列表

#### 返回值：

* 任意结构对象列表的Fasta格式字符串

#### 例：

``` Cpp
auto proPtr = new Protein;

auto fastaStr = fastaStr(proPtr->sub());
```

### 7.7 dumpFasta

``` Cpp
template <typename T>
void dumpFasta(const T &structPtrList, const string &dumpFilePath,
    const string &fileMode = "w")
```

将任意结构对象列表输出至Fasta文件。

#### 参数：

* structPtrList：任意结构对象列表
* dumpFilePath：输出文件路径
* fileMode：文件打开模式

#### 返回值：

* void

#### 例：

``` Cpp
auto proPtr = new Protein;

dumpFasta(proPtr->sub(), "xxx.fasta");
```

## 8. 常量

### 8.1 DIH

``` Cpp
enum class DIH;
```

枚举变量，表示主链二面角种类。DIH::PHI或DIH::L表示Phi，DIH::PSI或DIH::R表示Psi。

### 8.2 SIDE

``` Cpp
enum class SIDE;
```

枚举变量，表示主链二面角旋转时的转动侧。SIDE::N或SIDE::L表示转动N端，SIDE::C或SIDE::R表示转动C端。

### 8.3 RESIDUE_NAME_THREE_TO_ONE_MAP, RESIDUE_NAME_ONE_TO_THREE_MAP

``` Cpp
const unordered_map<string, string> RESIDUE_NAME_THREE_TO_ONE_MAP;
const unordered_map<string, string> RESIDUE_NAME_ONE_TO_THREE_MAP;
```

三字母，单字母残基名的相互转换哈希表。

## 9. 补充说明

### 9.1 解析函数

* PDB文件解析函数（load、loadModel）将完全按照PDB文件内容进行解析，不会对结构进行任何排序、合并或重组操作
* Load函数在解析时会跳过任何非"ATOM"关键词开头的行（包括"MODEL"）；而LoadModel函数会跳过任何非"ATOM"或"MODEL"关键词开头的行
* 解析时会去除所有字符串类型属性双端的空格字符

### 9.2 对于创建新对象的判定

#### 9.2.1 load函数：

* Protein：只会在解析开始前创建唯一的一个，并最终返回这个对象
* Chain：解析开始时，以及每次检测到链名发生变化时（从上一个"ATOM"行到当前行），都会创建一个新的链对象
* Residue：解析开始时，创建新链时，以及每次检测到残基名、残基编号或残基插入字符三者之一发生变化时（从上一个"ATOM"行到当前行），都会创建一个新的残基对象
* Atom：每检测到一个新的"ATOM"行都会创建一个新的Atom对象

#### 9.2.2 loadModel函数：

* Protein：解析开始前，以及每次检测到"MODEL"关键词时，都会创建一个新的蛋白对象。如果解析开始前创建的这个蛋白对象在函数返回前仍然为空，则其将在函数返回前被删除并析构
* Chain：解析开始时，创建新Model时，以及每次检测到链名发生变化时（从上一个"ATOM"行到当前行），都会创建一个新的链对象
* Residue：解析开始时，创建新Model时，创建新链时，以及每次检测到残基名、残基编号或残基插入字符三者之一发生变化时（从上一个"ATOM"行到当前行），都会创建一个新的残基对象
* Atom：每检测到一个新的"ATOM"行都会创建一个新的Atom对象
