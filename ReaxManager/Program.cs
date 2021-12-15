using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace ReactionManager2
{
    class ReactionManager
    {
        public List<List<List<int>>> atomList;
        public List<List<int>> a2mList;
        public List<List<List<int>>> m2aList;
        public List<int> idToType;
        public Dictionary<int, string> typeToAtom;
        public List<List<string>> molNumToString;
        public List<List<string>> molNumToSmiles;
        public ReactionManager()
        {
            atomList = new List<List<List<int>>>();
            a2mList = new List<List<int>>();
            m2aList = new List<List<List<int>>>();
            idToType = new List<int>();
            typeToAtom = new Dictionary<int, string>();
        }

        static void Main(string[] args)
        {
            ReactionManager reactionManager = new ReactionManager();
            reactionManager.FileOpen("bondsH2OO2_short.reaxc", 48268, 252);
            reactionManager.GetDataPerTime(0, 252);
        }


        public void A2mList(List<int> aList, List<int> a2m, List<List<int>> m2a, List<List<int>> mList, List<int> molNumList, List<int> xAList, List<int> progressList)
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

        public void SolveReaction(List<int> aList, List<int> a2mR, List<int> a2mP, List<List<int>> m2aR, List<List<int>> m2aP,
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

        public void FileOpen(string filePath, int atomNum, int timeStep)
        {
            atomList = new List<List<List<int>>>();
            Console.WriteLine(atomList.Count);
            for (int i = 0; i < timeStep; i++)
            {
                atomList.Add(new List<List<int>>());
                for (int j = 0; j < atomNum; j++)
                {
                    atomList[i].Add(new List<int>());
                }
            }
            idToType = new List<int>(atomNum);
            for (int i = 0; i < atomNum; i++)
            {
                idToType.Add(0);
            }
            typeToAtom = new Dictionary<int, string>()
            {   
                {1,"O"},{2,"O"},{3,"H"},{4,"C"},{5,"C"},{6,"O"},{7,"C"},{8,"O"},
                {9,"F"},{10,"H"},{11,"H"},{12,"_C_"},{13,"_C_"},{14,"_C_"},{15,"H"},{16,"_C_"},{17,"_C_"},{18,"_C_"},{19, "_C_"}, {20, "H"}, {21, "_C_"}
            };

            StreamReader streamReader = new StreamReader(filePath);
            StreamReader sr = streamReader;
            int timeStepCount = -1;
            while (sr.EndOfStream == false)
            {
                string[] line = sr.ReadLine().Split(" ");
                if (line[0] == "#")
                {
                    continue;
                }
                if (line[1] == "1")
                {                    
                    timeStepCount++;
                }
                //Console.WriteLine($"{atomList[timeStepCount].Count}");
                int atomID = int.Parse(line[1]);
                for (int j = 0; j < int.Parse(line[3]); j++)
                {
                    atomList[timeStepCount][atomID - 1].Add(int.Parse(line[4 + j]));
                }
                idToType[atomID - 1] = int.Parse(line[2]);
            }
        }

        public string ChangeFromIDToString(List<int> aList, int timeStep, int targetTime)
        {
            string atomChain = "";
            foreach (int ID in aList)
            {
                if (ID < 0)
                {
                    atomChain += "(";
                    continue;
                }
                if (ID >= 50000)
                {
                    atomChain += ")";
                    continue;
                }
                int type = idToType[ID - 1];
                string atomName = typeToAtom[type];
                if (timeStep == targetTime)
                {
                    if (!IsListEquall(this.atomList[timeStep][ID - 1], this.atomList[timeStep + 1][ID - 1]))
                    {
                        if (atomName == "O" && this.atomList[timeStep][ID - 1].Count == 1)
                        {
                            atomChain += "=";
                        }
                        atomChain += ("#" + atomName + "#");
                        continue;
                    }
                }
                else
                {
                    if (!IsListEquall(this.atomList[timeStep][ID - 1], this.atomList[targetTime][ID - 1]))
                    {
                        if (atomName == "O" && this.atomList[timeStep][ID - 1].Count == 1)
                        {
                            atomChain += "=";
                        }
                        atomChain += ("#" + atomName + "#");
                        continue;
                    }
                }
                if (atomName == "O" && this.atomList[timeStep][ID - 1].Count == 1)
                {
                    atomChain += "=";
                }
                if (atomName == "H")
                {
                    atomChain += "X";
                    continue;
                }
                atomChain += atomName;
            }
            return atomChain;
        }

        public bool IsListEquall(List<int> a, List<int> b)
        {
            if (a.Count != b.Count)
            {
                return false;
            }
            foreach (int aValue in a)
            {
                if (!b.Contains(aValue))
                {
                    return false;
                }
            }
            return true;
        }


        //ある原子IDを受け取って、その原子と結合している原子の中から、一番繋がっている原子IDを返す。(目的: 鎖の特定)
        public int FindMostChainedAtom(int atomID, List<List<int>> atomList_copy, int timeStep, bool isRemoved)
        {
            if (atomList_copy[atomID - 1] == null)
            {
                return -1;//結合している全ての原子を列挙済み
            }          
            int atom = atomList_copy[atomID - 1][0];
            foreach (int _atom in atomList_copy[atomID - 1])
            {
                if (atomList[timeStep][_atom - 1].Count > atomList[timeStep][atom - 1].Count)
                {
                    atom = _atom;
                }
            }
            if(!isRemoved)
            {
                return atom;
            }
            atomList_copy[atomID - 1].Remove(atom);
            atomList_copy[atom - 1].Remove(atomID);
            if (atomList_copy[atomID - 1].Count == 0)
            {
                atomList_copy[atomID - 1] = null;
            }
            if (atomList_copy[atom - 1].Count == 0)
            {
                atomList_copy[atom - 1] = null;
            }
            return atom;
        }

        public void SolveAllAtomsReaction(int time, bool HBIsIgnored, ref int totalMolNum, ref int totalSpecies, ref int totalD4OHNum, ref int totalH2ONum, ref int totalO2Num, ref int totalReactionNum, List<List<List<List<string>>>> totalReaction, List<Dictionary<string, int>> totalSpeciesDict)
        {
            molNumToString = new List<List<string>>();
            molNumToSmiles = new List<List<string>>();
            m2aList = new List<List<List<int>>>();
            for (int i = 0; i < 2; i++)
            {
                m2aList.Add(new List<List<int>>());
            }
            a2mList = new List<List<int>>();
            for (int i = 0; i < 2; i++)
            {
                a2mList.Add(new List<int>());
                for (int j = 0; j < atomList[0].Count; j++)
                {
                    a2mList[i].Add(0);
                }
            }
            for (int i = time; i <= time + 1; i++)
            {
                for (int j = 0; j < atomList[i].Count; j++)
                {
                    if ((typeToAtom[idToType[j]] == "H" || typeToAtom[idToType[j]] == "F") && atomList[i][j].Count >= 2 && HBIsIgnored)
                    {
                        int mostChainedID = FindMostChainedAtom(j + 1, atomList[i], i, false);                       
                        for (int k = 0; k < atomList[i][j].Count; k++)
                        {
                            if(atomList[i][j][k] == mostChainedID)
                            {
                                continue;
                            }
                            atomList[i][atomList[i][j][k] - 1].Remove(j + 1);
                        }
                        atomList[i][j] = new List<int>() { mostChainedID };
                    }
                    for (int k = 0; k < atomList[i][j].Count; k++)
                    {
                        if ((typeToAtom[idToType[j]] == "_C_" && typeToAtom[idToType[atomList[i][j][k] - 1]] != "_C_") || (typeToAtom[idToType[j]] != "_C_" && typeToAtom[idToType[atomList[i][j][k] - 1]] == "_C_"))
                        {
                            atomList[i][j][k] = -1;
                        }
                    }
                    atomList[i][j].RemoveAll(item => item == -1);
                }
            }

            for (int i = 0; i < 2; i++)
            {
                molNumToString.Add(new List<string>());
            }
            for (int i = 0; i < 2; i++)
            {
                molNumToSmiles.Add(new List<string>());
            }
            for (int i = 0; i < 2; i++)
            {
                int j = 0;
                List<List<int>> atomList_copy = new List<List<int>>(atomList[time + i]);
                for (int k = 0; k < atomList[time + i].Count; k++)
                {
                    atomList_copy[k] = new List<int>(atomList[time + i][k]);
                    if (atomList_copy[k].Count == 0 || typeToAtom[idToType[k]] == "_C_")
                    {
                        atomList_copy[k] = null;
                    }
                }
                while (true) //atomList_copyが全てnullになるまで繰り返す。(全ての分子の完成)
                {
                    if(atomList_copy[j] == null)
                    {
                        int j_count2 = 0;
                        for (int k = 0; k < atomList_copy.Count; k++)
                        {
                            j_count2++;
                            if (atomList_copy[k] != null)
                            {
                                j = k;
                                break;
                            }
                        }

                        if (j_count2 >= atomList_copy.Count)
                        {
                            break;
                        }
                    }
                    List<int> molList_temp = new List<int>() { j + 1 };
                    List<int> molList_temp2 = new List<int>() { j + 1 };
                    List<int> preMolList_temp = new List<int>();
                    int chainCount = 0;
                    while (molList_temp.Count != preMolList_temp.Count) // 追加できる原子idがなくなるまで繰り返す。(1分子の完成)
                    {
                        preMolList_temp = new List<int>(molList_temp);
                        foreach (int atom in preMolList_temp)//preMolListの各原子について、繰り返す。(側鎖の導出)
                        {
                            if (atom < 0 || atom >= 50000)
                            {
                                continue;
                            }
                            int atomID = atom;
                            List<int> chainedList = new List<int>();
                            while (true) //鎖の末端まで繰り返す。(一本の鎖の完成)
                            {
                                int addAtom = FindMostChainedAtom(atomID, atomList_copy, time + i, true); //チェインリストに加える原子を決める。(atomList_copyの結合情報も削除)
                                if (addAtom == -1)
                                {
                                    break;
                                }
                                else
                                {
                                    if (molList_temp.Contains(addAtom) || chainedList.Contains(addAtom))
                                    {
                                        break;
                                    }
                                    atomID = addAtom;
                                    chainedList.Add(addAtom);

                                }
                            }

                            if (chainCount == 0)
                            {
                                molList_temp.InsertRange(molList_temp.IndexOf(atom) + 1, chainedList);
                                chainedList.RemoveAll(item => typeToAtom[idToType[item-1]] == "H");
                                molList_temp2.InsertRange(molList_temp2.IndexOf(atom) + 1, chainedList);
                            }
                            else if (chainCount != 0 && chainedList.Count != 0)
                            {
                                molList_temp.Insert(molList_temp.IndexOf(atom) + 1, -1 * atom);
                                molList_temp.InsertRange(molList_temp.IndexOf(atom) + 2, chainedList);
                                molList_temp.Insert(molList_temp.IndexOf(atom) + 2 + chainedList.Count, 50000 + atom);
                                if (chainedList.Count == 1 && typeToAtom[idToType[chainedList[0] - 1]] == "H")
                                {
                                    continue;
                                }
                                molList_temp2.Insert(molList_temp2.IndexOf(atom) + 1, -1 * atom);
                                molList_temp2.InsertRange(molList_temp2.IndexOf(atom) + 2, chainedList);
                                molList_temp2.Insert(molList_temp2.IndexOf(atom) + 2 + chainedList.Count, 50000 + atom);
                            }
                        }
                        chainCount++;
                    }
                    string smiles_temp = ChangeFromIDToString(molList_temp2, time + i, time).Replace("X", "");
                    if (smiles_temp != "" && smiles_temp != "#H#" )
                    {
                        molNumToSmiles[i].Add(smiles_temp);
                        m2aList[i].Add(molList_temp2);
                    }                  
                    if (molList_temp.Count != 0)
                    {
                        molNumToString[i].Add(ChangeFromIDToString(molList_temp, time + i, time));                                                                        
                    }
                    
                }
            }
            totalMolNum = molNumToSmiles[0].Count;
            List<string> Smiles = new List<string>(molNumToSmiles[0]);
            List<int> smilesCount = new();
            Dictionary<string, int> stringToSmilesCount = new();
            for (int i=0; i<Smiles.Count; i++)
            {
                Smiles[i] = Smiles[i].Replace("#", "");
                Smiles[i] = Smiles[i].Replace("H", "");
                string Smiles2 = Smiles[i];
                Smiles[i] = Smiles[i].Replace("=", "");
                Smiles[i] = Smiles[i].Replace("(", "");
                Smiles[i] = Smiles[i].Replace(")", "");
                if(!smilesCount.Contains(Smiles[i].Count()) && Smiles[i] != "=O")
                {
                    smilesCount.Add(Smiles[i].Count());
                    stringToSmilesCount.Add(Smiles2,1);
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

            //それぞれの分子を特定して、数え上げていく
            List<string> molNumToSmiles_copy = new(molNumToSmiles[0]);
            foreach (string smiles in molNumToSmiles_copy)
            {
                if (smiles.Length > 1000)
                {
                    continue;
                }

                string smiles2 = smiles.Replace("H", "");
                smiles2 = smiles2.Replace("#", "");
                smiles2 = smiles2.Replace("(", "");
                smiles2 = smiles2.Replace(")", "");
                smiles2 = smiles2.Replace("=", "");
                switch (smiles2.Length)
                {
                    case 1:
                        totalH2ONum += 1;
                        break;
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
            for (int i = 0; i < totalTime - 1; i++)
            {
                int molNum = 0;
                int species = 0;
                int D4OHNum = 0;
                int H2ONum = 0;
                int O2Num = 0;
                int ReactionNum = 0;
                SolveAllAtomsReaction(i, true, ref molNum,ref species, ref D4OHNum, ref H2ONum, ref O2Num, ref ReactionNum, totalReaction, totalSpeciesDict);
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

            File.AppendAllText(@"ReactionData.txt", "List of Species" + Environment.NewLine);
            for (int i = 0; i < totalTime - 1; i += totalTime - 2)
            {
                File.AppendAllText(@"ReactionData.txt", $"TimeStep{i}" + Environment.NewLine);
                foreach (KeyValuePair<string, int> smiles in totalSpeciesDict[i])
                {
                    File.AppendAllText(@"ReactionData.txt", $"●{smiles.Key}: {smiles.Value}" + Environment.NewLine);
                }
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
