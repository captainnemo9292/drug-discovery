# drug-discovery

Graph Neural Network based Generative Modeling for COVID-19 protease inhibitor Drug Discovery 
- by 류병우

모바일 어플리케이션: (직접 해보세요!!!) ./drug-discovery.apk

데모 영상: 

모델 & 데이터셋: 

안녕하세요 Github, 딥러닝에 관심있는 고등학생입니다. 전세계적으로 코로나바이러스가 확산됨에 따라 AI 커뮤니티에서도 컴퓨터 과학을 기반으로 이 위기를 해결하고자 하는 연구 활동이 활발하게 진행되고 있습니다. 이러한 코로나 퇴치를 위한 노력에 영감을 받아 저 또한 그래프 인공신경망과 Generative Modeling을 신약 후보 물질 개발에 적용하는 프로젝트를 실시하게 되었습니다. 다음은 GNN 기반 모델들을 기반으로 코로나 증식 효소 억제 물질 생성 및 ligand와 protein의 결합 가능성 예측에 도움을 주는 모바일 어플리케이션을 개발하는 토이프로젝트의 과정입니다. 한 번씩 봐주시면 감사하겠습니다.

1. 신약 후보 물질 데이터 확보 

신약 후보 물질 구조 생성을 위해서는 인공지능을 훈련시키기 전에 코로나 바이러스 복제 효소를 억제하는 특성을 지닌 화학 물질 학습 데이터가 필요했습니다. 따라서, PDB (Protein Data Bank: 분자의 3D 구조를 제공하는 공공 데이터베이스)에서 바이러스의 증식을 촉진하는 효소를 검색하여 단백질 분해 효소 6LU7이 COVID-19의 증식을촉진시킨다는사실을 알아냈습니다. 이후 virtual screening 플랫폼인 Pharmit 에서 6LU7 효소와 결합 가능한 물질을 알아보았습니다. Pharmit에서 Docking 시뮬레이션을 가동한 결과, binding affinity가 특히 월등한 화학 물질 약 25000 개의 SMILES 데이터를 ZINC 화학 분자 데이터베이스에서 추출할 수 있었습니다.

2. 그래프 인공신경망을 활용한 모델 학습

화학 구조를 보다 더 효과적으로 추상화하기 위해 SMILES 데이터를 그래프로 변환하여 (원자를 node, 화학 결합을 edge로 표현한 네트워크 데이터) 그레프 인공신경망의 일종이자 Generative 모델에 해당하는 DGMG (Deep Generative Models of Graphs) 모델에 학습시켜보았습니다. 또한, PDBBind 데이터셋을 ACNN (Atomic Convolutional Networks) 모델에 학습시켜 ligand 분자와 protein 분자가 주어지면 둘의 3D 구조를 고려하여 binding affinity를 예측하는 모델을 개발했습니다. 

3. 모바일 어플리케이션 개발

학습시킨 모델들은 Flask 웹프레임워크 기반의 API로 배포하였고, 그러한 API를 활용하여 신약 후보 물질 개발과 virtual screening 서비스를 제공하는 모바일 어플리케이션을 개발해 보았습니다. React Native와 Expo를 활용하여 구축한 앱은 DGMG 모델이 생성한 분자와 유사한 신약 후보 물질들의 정보를 PubChem 데이터베이스에서 찾아 풍부한 메타데이터와 함께 반환해줍니다. (메타데이터: SMILES 데이터, 각종 정식 명칭 및 화학 수치 자료 등 PubChem 데이터베이스 API에서 제공하는 화학적 특성 정보) 또한, 생성된 ligand 분자와 사용자가 입력한 protein 분자 간의 binding affinity를 ACNN 모델로 계산하고, toxicity 여부 또한 Deep Graph Library 개발진이 제공하는 GCN 모델로 예측하여 virtual screening 분석 결과를 제공해줍니다.  

p.s. 과학적인 규범에 따라 진행되지 않았고 전문가가 아닌 개인에 의해 이루어진 토이프로젝트임을 가만해 주십시오. 프로젝트를 알리는 이유는 관심사가 저와 비슷한 분들께 프로젝트 아이디어를 갖게 하게끔 하거나, 저같은 학생 분들이 Generative Modeling 과 GNN이 신약 개발에도 적용될 수 있다는 사실을 알게끔 하고 싶어서 입니다.
