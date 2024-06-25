# <p align="center">![image](./png/ico.png) Dear-DIA<sup>XMBD</sup></p> 
## <p align="center">Deep autoencoder for data-independent acquisition proteomics</p> 

## Introduction
Dear-DIA<sup>XMBD</sup> is a spectrum-centric method that combines the deep variational autoencoder (VAE) with other machine learning algorithms, to detect the correspondence between precursors and fragments in DIA data without the help of DDA experiments. 

- Dear-DIA<sup>XMBD</sup> produces the pseudo-tandem MS spectra to search database, and generates the internal libraries. 
- Dear-DIA<sup>XMBD</sup> can be easily integrated into the existing workflow, because the output file of Dear-DIA<sup>XMBD</sup> is in MGF format that can be processed by common search engines, including Comet, X!Tandem and MSFragger etc. 
- Furthermore, benefiting from the fact that the autoencoder is an unsupervised deep learning model, Dear-DIA<sup>XMBD</sup> shows an excellent performance on the DIA data of different species obtained by different instrument platforms.
- Dear-DIA<sup>XMBD</sup> is a **cross-platform open source** software and follows the GPL v3.0 licence. **Users can download and use Dear-DIA<sup>XMBD</sup> for free.**

**Download:** [Dear-DIA-GUI-win-x64-free-install.zip](https://github.com/jianweishuai/Dear-DIA-XMBD/releases)

## Installation
If you have any problems installing the software, please feel free to contact us.

### Windows System
**Free installation:** [Dear-DIA-GUI-win-x64-free-install.zip](https://github.com/jianweishuai/Dear-DIA-XMBD/releases)
```
1. download the newest "Dear-DIA-GUI-win-x64-free-install.zip" file.
2. decompress this file to any path.
3. double click "Dear-DIA-GUI.exe" with a blue butterfly icon..
```

### Linux System
**Linux source code:** [Dear-DIA-linux-src-x86_64.zip](https://github.com/jianweishuai/Dear-DIA-XMBD/releases)
- **Note: need g++ >= 5.4.0 or compiler support c++11**
- "Dear-DIA-linux-src-x64.zip" file contains source code, dependencies and binary release file “deardia” compiled on CentOS 7.8 system. Users can try to run the executable file “deardia” on other linux systems such as Ubuntu system.
```
1. download the "Dear-DIA-linux-src-x86_64.zip" file.
2. decompress the "Dear-DIA-linux-src-x86_64.zip" to your linux server.
3. Enter the "make" command in the unzipped directory containing the Makefile.
4. Complete compilation.
```

## Usage

### Windows 
```
1. Double click "Dear-DIA-GUI.exe"
```
![image](./png/1.PNG)

```
2. Configure the parameters
```
![image](./png/2.PNG)

```
3. Click "start" button and wait for finishing.
```
![image](./png/3.PNG)

### Linux

```
./deardia --config=deardia.config.new --out_dir=/path/to/ouput_dir --input=/path/to/mzXML_file
```

## Update News

- **2022.6.01 updatde:** Upgrade Dear-DIA to support mzML, mzXML, wiff and raw formats.
- **2022.5.30 updatde:** fixd some BUGs.

## Supported data format

- **2022.6.01 updatde:** Support wiff (Sciex) and raw (ThermoFisher) formats.
- **2022.5.31 updatde:** Support mzML format.
- **2022.5.21 updatde:** Support mzXML format.

## Future updates
We will consider the following features in future releases.
- Support vendor formats including .wiff and .raw.
- Support more post-translation modifications.

## Contact

- Please post any questions, feedback, comments or suggestions on the  [GitHub Discussion board](https://github.com/jianweishuai/Dear-DIA-XMBD/issues).
- Email: qingzuhe@stu.xmu.edu.cn; jianweishuai@xmu.edu.cn; jhan@xmu.edu.cn
