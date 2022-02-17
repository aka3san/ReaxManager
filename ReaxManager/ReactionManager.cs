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
        public ReactionManager(AtomInputData atomInputData, MoleculeManager moleculeManager)
        {
            this.atomInputData = atomInputData;            
        }

        static void Main(string[] args)
        { 
            /*
            ReactionManager reactionManager = new ReactionManager();
            reactionManager.FileOpen("bondsH2OO2_short.reaxc", 48268, 252);
            reactionManager.GetDataPerTime(0, 252);
            */
            
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
                if(typeToAtom[idToType[aList[i]-1]] == "H")
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

        void SolveReaction(List<int> aList, List<int> a2mR, List<int> a2mP, List<List<int>> m2aR, List<List<int>> m2aP,
                                    List<List<int>> reacMList, List<List<int>> prodMList, List<int> molListReac, List<int> molListProd, int progress, List<int> reacProgressList, List<int> prodProgressList, List<bool> isChanged, int time)
        {
            if (IsListEquall(atomList[time][progress], atomList[time + 1][progress]) || typeToAtom[idToType[progress]] == "H" || typeToAtom[idToType[progress]] == "O")
            {
                return;
            }           
            else
            {
                isChanged[progress] = true;
            }

            if (reacProgressList.Contains(a2mR[progress]) || prodProgressList.Contains(a2mP[progress]))
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
                A2mList(list1, a2mP, m2aP, prodMList, molListProd, list2, prodProgressList);
                A2mList(list2, a2mR, m2aR, reacMList, molListReac, list3, reacProgressList);
                list2.Clear();
            }
            //Console.WriteLine($"進行度: {progress}/{ 48168}");
        }

        

        

        

        

        public void SolveAllAtomsReaction(int time, bool HBIsIgnored, ref int totalMolNum, ref int totalSpecies, ref int totalD4OHNum, ref int totalH2ONum, ref int totalO2Num, ref int totalReactionNum, List<List<List<List<string>>>> totalReaction, List<Dictionary<string, int>> totalSpeciesDict, List<Dictionary<string,string>> totalSmilesPlusH)
        {

            MoleculeManager moleculeManager = new MoleculeManager(atomInputData, time);
            totalMolNum = moleculeManager.MolNumToSmiles[0].Count;
            List<string> Smiles = new List<string>(moleculeManager.MolNumToSmiles[0]);
            List<int> smilesCount = new();
            Dictionary<string, int> stringToSmilesCount = new();
            Dictionary<string, string> smilesPlusH = new();
            for (int i=0; i<Smiles.Count; i++)
            {
                Smiles[i] = Smiles[i].Replace("[", "");
                Smiles[i] = Smiles[i].Replace("]", "");
                Smiles[i] = Smiles[i].Replace(":", "");
                Smiles[i] = Smiles[i].Replace("1", "");
                Smiles[i] = Smiles[i].Replace("H", "");
                if(Smiles[i] == "=O")
                {
                    continue;
                }
                string Smiles2 = Smiles[i];
                Smiles[i] = Smiles[i].Replace("=", "");
                Smiles[i] = Smiles[i].Replace("(", "");
                Smiles[i] = Smiles[i].Replace(")", "");
                if(!smilesCount.Contains(Smiles[i].Count()))
                {
                    smilesCount.Add(Smiles[i].Count());
                    stringToSmilesCount.Add(Smiles2,1);
                    smilesPlusH.Add(Smiles2, molNumToString[0][i]);
                }
                else
                {
                    foreach(string key in stringToSmilesCount.Keys)
                    {
                        if(Smiles[i].Length == key.Replace("=","").Replace("(","").Replace(")","").Length)
                        {
                            stringToSmilesCount[key] += 1;
                        }
                    }
                }
            }
            totalSpecies = smilesCount.Count;
            totalSpeciesDict.Add(stringToSmilesCount);
            totalSmilesPlusH.Add(smilesPlusH);

            //それぞれの分子を特定して、数え上げていく
            List<string> molNumToSmiles_copy = new(molNumToSmiles[0]);
            foreach (string smiles in molNumToSmiles_copy)
            {
                if (smiles.Length > 1000)
                {
                    continue;
                }

                string smiles2 = smiles.Replace("H", "");
                smiles2 = smiles2.Replace("[", "");
                smiles2 = smiles2.Replace("]", "");
                smiles2 = smiles2.Replace(":", "");
                smiles2 = smiles2.Replace("1", "");
                smiles2 = smiles2.Replace("(", "");
                smiles2 = smiles2.Replace(")", "");
                if (smiles2.Length == 1)
                {
                    totalH2ONum += 1;
                }
                smiles2 = smiles2.Replace("=", "");
                switch (smiles2.Length)
                {                   
                    case 2:
                        totalO2Num += 1;
                        break;
                    case 111:
                        totalD4OHNum += 1;
                        break;
                    case 113:
                        totalD4OHNum += 1;
                        break;

                }
                if (smiles2.Length % 111 == 0 && smiles2.Length > 111)
                {
                    totalD4OHNum += smiles2.Length / 111 - 1;
                }
            }

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < m2aList[i].Count; j++)
                {
                    foreach (int a in m2aList[i][j])
                    {
                        if (a < 0 || a >= 50000)
                        {
                            continue;
                        }
                        a2mList[i][a - 1] = j + 1;
                    }
                }
            }
            List<bool> IsChanged = new();
            for (int i = 0; i < atomList[0].Count; i++)
            {
                IsChanged.Add(false);
            }
            List<int> reacProgressList = new List<int>();
            List<int> prodProgressList = new List<int>();
            List<List<List<List<int>>>> reac_prod_List = new List<List<List<List<int>>>>(atomList[0].Count);
            for (int i = 0; i < atomList[0].Count; i++)
            {
                reac_prod_List.Add(new List<List<List<int>>>());
                for (int j = 0; j < 2; j++)
                {
                    reac_prod_List[i].Add(new List<List<int>>());
                }
            }
            List<List<List<List<string>>>> reac_prod_List_str = new List<List<List<List<string>>>>(atomList[0].Count);
            for (int i = 0; i < atomList[0].Count; i++)
            {
                reac_prod_List_str.Add(new List<List<List<string>>>());
                for (int j = 0; j < 2; j++)
                {
                    reac_prod_List_str[i].Add(new List<List<string>>());
                }
            }
            List<List<List<string>>> reac_prod_StringList = new();
            for (int i = 0; i < atomList[0].Count; i++)
            {
                reac_prod_StringList.Add(new List<List<string>>());
                for (int j = 0; j < 2; j++)
                {
                    reac_prod_StringList[i].Add(new List<string>());
                }
            }
            List<List<int>> reacMList = new List<List<int>>(atomList[0].Count);
            for (int i = 0; i < atomList[0].Count; i++)
            {
                reacMList.Add(new List<int>());
            }
            List<List<int>> prodMList = new List<List<int>>(atomList[0].Count);
            for (int i = 0; i < atomList[0].Count; i++)
            {
                prodMList.Add(new List<int>());
            }
            List<List<int>> atomList2 = new List<List<int>>(atomList[time]);

            for (int i = 0; i < atomList[time].Count; i++)
            {
                SolveReaction(atomList[time][i], a2mList[0], a2mList[1], m2aList[0], m2aList[1],
                              reac_prod_List[i][0], reac_prod_List[i][1], reacMList[i], prodMList[i], i, reacProgressList, prodProgressList, IsChanged, time);
                if (reacMList[i].Count != 0 && IsChanged[i])
                {
                    foreach (int molNum in reacMList[i])
                    {
                        if (molNum == 0)
                        {
                            continue;
                        }
                        string molStr = molNumToSmiles[0][molNum - 1];
                        reac_prod_StringList[i][0].Add(molStr);
                    }
                    foreach (int molNum in prodMList[i])
                    {
                        if (molNum == 0)
                        {
                            continue;
                        }
                        string molStr = molNumToSmiles[1][molNum - 1];
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
            Console.WriteLine($"reac_prod_StringList.Count: {reac_prod_StringList.Count}");
        }

        public void GetDataPerTime(int molID, int totalTime)
        {
            List<int> totalMolNum = new();
            List<int> totalSpecies = new();
            List<int> totalD4OHNum = new();
            List<int> totalH2ONum = new();
            List<int> totalO2Num = new();
            List<int> totalReactionNum = new();
            List<List<List<List<string>>>> totalReaction = new();
            List<Dictionary<string, int>> totalSpeciesDict = new();
            List<Dictionary<string, string>> totalSmilesPlusH = new();
            for (int i = 0; i < totalTime - 1; i++)
            {
                int molNum = 0;
                int species = 0;
                int D4OHNum = 0;
                int H2ONum = 0;
                int O2Num = 0;
                int ReactionNum = 0;
                SolveAllAtomsReaction(i, true, ref molNum,ref species, ref D4OHNum, ref H2ONum, ref O2Num, ref ReactionNum, totalReaction, totalSpeciesDict, totalSmilesPlusH);
                totalMolNum.Add(molNum);
                totalSpecies.Add(species);
                totalD4OHNum.Add(D4OHNum);
                totalH2ONum.Add(H2ONum);
                totalO2Num.Add(O2Num);
                totalReactionNum.Add(ReactionNum);
                Console.WriteLine($"進行度: {i + 1}/{totalTime}");
            }            
            Console.WriteLine("テキストファイルが出力されました。");

            //テキストファイルに出力
            try
            {
                File.WriteAllText(@"ReactionData.txt", "Reaction Data" + Environment.NewLine);
            }
            catch (IOException e)
            {
                Console.WriteLine(e.Message);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of Molecules" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalMolNum[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of D4OH" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalD4OHNum[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of H2O" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalH2ONum[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of O2" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalO2Num[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of Reaction" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalReactionNum[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "Number of Species" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{totalSpecies[i]}" + Environment.NewLine);
            }


            List<int> polymolCountPerTime = new List<int>();
            File.AppendAllText(@"ReactionData.txt", "List of Species" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                //重合した分子の数
                int polymolCount = 0;

                File.AppendAllText(@"ReactionData.txt", $"TimeStep{i}" + Environment.NewLine);
                foreach (KeyValuePair<string, int> smiles in totalSpeciesDict[i])
                {
                    string molFormula = "";
                    int cCount, hCount, fCount, oCount;
                    cCount = 0; hCount = 0; fCount = 0; oCount = 0;
                    foreach(char atomName in totalSmilesPlusH[i][smiles.Key])
                    {
                        switch(atomName)
                        {
                            case 'C':
                                cCount++;
                                break;
                            case 'H':
                                hCount++;
                                break;
                            case 'F':
                                fCount++;
                                break;
                            case 'O':
                                oCount++;
                                break;
                        }
                    }
                    if(cCount > 37)
                    {
                        polymolCount++;
                    }
                    molFormula += (cCount != 0) ? ("C" + cCount.ToString()) : "";
                    molFormula += (hCount != 0) ? ("H" + hCount.ToString()) : "";
                    molFormula += (fCount != 0) ? ("F" + fCount.ToString()) : "";
                    molFormula += (oCount != 0) ? ("O" + oCount.ToString()) : "";
                    if(i == 0 || i== totalTime-2)
                    {
                        File.AppendAllText(@"ReactionData.txt", $"●{smiles.Key} \"({molFormula})\": {smiles.Value}" + Environment.NewLine);
                    }                   
                }
                polymolCountPerTime.Add(polymolCount);             
            }

            File.AppendAllText(@"ReactionData.txt", $"●重合した分子の数" + Environment.NewLine);
            for (int i=0; i<totalTime-1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"{polymolCountPerTime[i]}" + Environment.NewLine);
            }

            File.AppendAllText(@"ReactionData.txt", "ReactionPerTime" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i++)
            {
                File.AppendAllText(@"ReactionData.txt", $"TimeStep{i}" + Environment.NewLine);
                for(int i2=0; i2<totalReaction[i].Count; i2++)
                {
                    File.AppendAllText(@"ReactionData.txt", $"------反応{i2+1}------" + Environment.NewLine);
                    for(int j=0; j<totalReaction[i][i2].Count; j++)
                    {
                        for(int k=0; k<totalReaction[i][i2][j].Count; k++)
                        {
                            File.AppendAllText(@"ReactionData.txt", $"{totalReaction[i][i2][j][k]}" + Environment.NewLine);
                            if(k == totalReaction[i][i2][j].Count-1)
                            {
                                break;
                            }
                            File.AppendAllText(@"ReactionData.txt", "+" + Environment.NewLine);
                        }
                        if(j==0)
                        {
                            File.AppendAllText(@"ReactionData.txt", "--↓↓--" + Environment.NewLine);
                        }                       
                    }
                    File.AppendAllText(@"ReactionData.txt", "------------" + Environment.NewLine);
                }
            }
        }
    }
}
