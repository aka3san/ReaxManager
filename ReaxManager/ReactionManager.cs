using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace ReaxManager
{
    class ReactionManager
    {
        private AtomInputData atomInputData;       
        public ReactionManager(AtomInputData atomInputData)
        {
            this.atomInputData = atomInputData;            
        }        

        void A2mList(List<int> aList, List<int> a2m, List<List<int>> m2a, List<List<int>> mList, List<int> molNumList, List<int> xAList, List<int> progressList)
        {
            mList.Clear();
            xAList.Clear();
            molNumList.Clear();
            HashSet<int> molnum = new HashSet<int>();
            for (int i = 0; i < aList.Count; i++)
            {
                if (aList[i] < 0 || aList[i] >= 50000)
                {
                    continue;
                }
                if(atomInputData.TypeToAtom[atomInputData.IdToType[aList[i]-1]] == "H")
                {
                    continue;
                }
                molnum.Add(a2m[aList[i] - 1]);
                progressList.Add(a2m[aList[i] - 1]);
            }
            molNumList.AddRange(molnum);
            M2aList(molnum, m2a, mList);
            for (int i = 0; i < mList.Count; i++)
            {
                for (int j = 0; j < mList[i].Count; j++)
                {
                    xAList.Add(mList[i][j]);
                }
            }
        }

        public void M2aList(HashSet<int> molnum, List<List<int>> m2a, List<List<int>> mList)
        {
            foreach (int mol in molnum)
            {
                if (mol == 0)
                {
                    continue;
                }
                mList.Add(m2a[mol - 1]);
            }
        }

        void SolveReaction(List<int> aList, MoleculeManager moleculeManager, List<List<int>> reacMList, List<List<int>> prodMList, List<int> molListReac, List<int> molListProd,
                           int progress, List<int> reacProgressList, List<int> prodProgressList, List<bool> isChanged, int time)
        {
            if (MoleculeManager.IsListEquall(atomInputData.AtomList[time][progress], atomInputData.AtomList[time + 1][progress]) || atomInputData.TypeToAtom[atomInputData.IdToType[progress]] == "H" || atomInputData.TypeToAtom[atomInputData.IdToType[progress]] == "O")
            {
                return;
            }           
            else
            {
                isChanged[progress] = true;
            }

            if (reacProgressList.Contains(moleculeManager.A2mList[0][progress]) || prodProgressList.Contains(moleculeManager.A2mList[1][progress]))
            {
                reacMList.Add(new List<int>() { 0 });
                prodMList.Add(new List<int>() { 0 });
                return;
            }

            reacMList.Add(aList);
            List<int> list1 = new List<int>();
            List<int> list2 = new List<int>();
            List<int> list3 = new List<int>();
            list3 = new List<int>(aList);

            while (list1.Count != list3.Count)
            {
                list1 = new List<int>(list3);
                A2mList(list1, moleculeManager.A2mList[1], moleculeManager.M2aList[1], prodMList, molListProd, list2, prodProgressList);
                A2mList(list2, moleculeManager.A2mList[0], moleculeManager.M2aList[0], reacMList, molListReac, list3, reacProgressList);
                list2.Clear();
            }
            //Console.WriteLine($"進行度: {progress}/{ 48168}");
        }
     
        public void SolveAllAtomsReaction(int time, bool HBIsIgnored, ref int totalMolNum, ref int totalSpecies, List<List<int>> totalTargetMoleculeNum, List<string> targetMoleculeSmilesList, ref int totalReactionNum, List<List<List<List<string>>>> totalReaction, List<Dictionary<string, int>> totalSpeciesDict, List<Dictionary<string,string>> totalSmilesPlusH)
        {

            MoleculeManager moleculeManager = new MoleculeManager(atomInputData, time);
            totalMolNum = moleculeManager.MolNumToSmiles[0].Count;

            SortSmiles sortSmiles = new SortSmiles(moleculeManager.MolNumToSmiles[0], moleculeManager.MolNumToString[0]);
            Dictionary<string, int> stringToSpeciesNum = new();
            Dictionary<string, string> smilesPlusH = new();
            sortSmiles.SortSpecies(ref totalSpecies, stringToSpeciesNum, smilesPlusH);            
            totalSpeciesDict.Add(stringToSpeciesNum);
            totalSmilesPlusH.Add(smilesPlusH);

            List<int> targetMoleculerNum = sortSmiles.CountTargetMoleculeNumbers(targetMoleculeSmilesList);
            totalTargetMoleculeNum.Add(targetMoleculerNum);
            

            //反応解析↓
            List<bool> IsChanged = new();
            for (int i = 0; i < atomInputData.TotalAtomNum; i++)
            {
                IsChanged.Add(false);
            }
            List<int> reacProgressList = new List<int>();
            List<int> prodProgressList = new List<int>();
            List<List<List<List<int>>>> reac_prod_List = new List<List<List<List<int>>>>(atomInputData.TotalAtomNum);
            for (int i = 0; i < atomInputData.TotalAtomNum; i++)
            {
                reac_prod_List.Add(new List<List<List<int>>>());
                for (int j = 0; j < 2; j++)
                {
                    reac_prod_List[i].Add(new List<List<int>>());
                }
            }
            List<List<List<List<string>>>> reac_prod_List_str = new List<List<List<List<string>>>>(atomInputData.TotalAtomNum);
            for (int i = 0; i < atomInputData.TotalAtomNum; i++)
            {
                reac_prod_List_str.Add(new List<List<List<string>>>());
                for (int j = 0; j < 2; j++)
                {
                    reac_prod_List_str[i].Add(new List<List<string>>());
                }
            }
            List<List<List<string>>> reac_prod_StringList = new();
            for (int i = 0; i < atomInputData.TotalAtomNum; i++)
            {
                reac_prod_StringList.Add(new List<List<string>>());
                for (int j = 0; j < 2; j++)
                {
                    reac_prod_StringList[i].Add(new List<string>());
                }
            }
            List<List<int>> reacMList = new List<List<int>>(atomInputData.TotalAtomNum);
            for (int i = 0; i < atomInputData.TotalAtomNum; i++)
            {
                reacMList.Add(new List<int>());
            }
            List<List<int>> prodMList = new List<List<int>>(atomInputData.TotalAtomNum);
            for (int i = 0; i < atomInputData.TotalAtomNum; i++)
            {
                prodMList.Add(new List<int>());
            }
            List<List<int>> atomList2 = new List<List<int>>(atomInputData.TotalAtomNum);

            for (int i = 0; i < atomInputData.TotalAtomNum; i++)
            {
                SolveReaction(atomInputData.AtomList[time][i], moleculeManager, reac_prod_List[i][0], reac_prod_List[i][1], reacMList[i], prodMList[i],
                              i, reacProgressList, prodProgressList, IsChanged, time);
                if (reacMList[i].Count != 0 && IsChanged[i])
                {
                    foreach (int molNum in reacMList[i])
                    {
                        if (molNum == 0)
                        {
                            continue;
                        }
                        string molStr = moleculeManager.MolNumToSmiles[0][molNum - 1];
                        reac_prod_StringList[i][0].Add(molStr);
                    }
                    foreach (int molNum in prodMList[i])
                    {
                        if (molNum == 0)
                        {
                            continue;
                        }
                        string molStr = moleculeManager.MolNumToSmiles[1][molNum - 1];
                        reac_prod_StringList[i][1].Add(molStr);
                    }
                }
                else
                {
                    reac_prod_StringList[i][0].Add("None");
                }
            }
            reac_prod_StringList.RemoveAll(item => item[0][0] == "None");
            totalReactionNum = reac_prod_StringList.Count;
            totalReaction.Add(reac_prod_StringList);            
        }

        public void GetDataPerTime()
        {
            List<int> totalMolNum = new();
            List<int> totalSpecies = new();
            List<List<int>> totalTargetMoleculeNum = new();      
            List<List<List<List<string>>>> totalReaction = new();
            List<Dictionary<string, int>> totalSpeciesDict = new();
            List<Dictionary<string, string>> totalSmilesPlusH = new();
            for (int i = 0; i < atomInputData.TotalTimeStep - 1; i++)
            {
                int molNum = 0;
                int species = 0;                
                int ReactionNum = 0;
                SolveAllAtomsReaction(i, true, ref molNum,ref species, totalTargetMoleculeNum, atomInputData.TargetMoleculeSmilesList, ref ReactionNum, totalReaction, totalSpeciesDict, totalSmilesPlusH);
                totalMolNum.Add(molNum);
                totalSpecies.Add(species);               
                totalReactionNum.Add(ReactionNum);
                Console.WriteLine($"進行度: {i + 1}/{atomInputData.TotalTimeStep}");
            }
            Console.WriteLine("テキストファイルが出力されました。");

            //テキストファイルに出力
            OutputReactionData outputReactionData = new OutputReactionData()
        }
    }
}
