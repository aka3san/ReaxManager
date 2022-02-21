using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ReaxManager
{
    public class MoleculeManager
    {
        private List<List<int>> a2mList;
        public List<List<int>> A2mList => a2mList;
        private List<List<List<int>>> m2aList;
        public List<List<List<int>>> M2aList => m2aList;
        private  List<List<string>> molNumToSmiles;
        public List<List<string>> MolNumToSmiles => molNumToSmiles;
        private List<List<string>> molNumToString;
        public List<List<string>> MolNumToString => molNumToString;

        private AtomInputData atomInputData;

        public MoleculeManager(AtomInputData atomInputData, int time)//ある時刻の分子一覧を管理するクラス
        {           
            this.atomInputData = atomInputData;
            molNumToSmiles = new List<List<string>>();
            for (int i = 0; i <2; i++)
            {
                molNumToSmiles.Add(new List<string>());
            }
            molNumToString = new List<List<string>>();
            for (int i = 0; i < 2; i++)
            {
                molNumToString.Add(new List<string>());
            }
            m2aList = new List<List<List<int>>>();
            for (int i = 0; i < 2; i++)
            {
                m2aList.Add(new List<List<int>>());
            }
            a2mList = new List<List<int>>();
            for (int i = 0; i < 2; i++)
            {
                a2mList.Add(new List<int>());
                for (int j = 0; j < atomInputData.AtomList[0].Count; j++)
                {
                    a2mList[i].Add(0);
                }
            }

            for (int i = 0; i < 2; i++)
            {
                int j = 0;
                List<List<int>> atomList_copy = new List<List<int>>(atomInputData.AtomList[time + i]);
                for (int k = 0; k < atomInputData.AtomList[time + i].Count; k++)
                {
                    atomList_copy[k] = new List<int>(atomInputData.AtomList[time + i][k]);
                    if (atomList_copy[k].Count == 0 || atomInputData.TypeToAtom[atomInputData.IdToType[k]] == "_C_")
                    {
                        atomList_copy[k] = null;
                    }
                }
                while (true) //atomList_copyが全てnullになるまで繰り返す。(全ての分子の完成)
                {
                    if (atomList_copy[j] == null)
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
                                int addAtom = atomInputData.RemoveMostChainedAtom(time + i, atomID, atomList_copy);                                
                                //チェインリストに加える原子を決める。(atomList_copyの結合情報も削除)
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
                                chainedList.RemoveAll(item => atomInputData.TypeToAtom[atomInputData.IdToType[item - 1]] == "H");
                                molList_temp2.InsertRange(molList_temp2.IndexOf(atom) + 1, chainedList);
                            }
                            else if (chainCount != 0 && chainedList.Count != 0)
                            {

                                molList_temp.Insert(molList_temp.IndexOf(atom) + 1, -1 * atom);
                                molList_temp.InsertRange(molList_temp.IndexOf(atom) + 2, chainedList);
                                molList_temp.Insert(molList_temp.IndexOf(atom) + 2 + chainedList.Count, 50000 + atom);
                                if (chainedList.Count == 1 && atomInputData.TypeToAtom[atomInputData.IdToType[chainedList[0] - 1]] == "H")
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
                    string smiles_temp = ChangeFromIDToString(molList_temp2, time + i, time, true).Replace("X", "");
                    string smiles_temp2 = ChangeFromIDToString(molList_temp, time + i, time, false);
                    if (smiles_temp != "" && smiles_temp != "[H:1]")
                    {
                        molNumToSmiles[i].Add(smiles_temp);
                        m2aList[i].Add(molList_temp2);
                    }
                    if (smiles_temp2 != "H" && smiles_temp2 != "[H:1]")
                    {
                        molNumToString[i].Add(smiles_temp2);
                    }                    
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
        }

        public string ChangeFromIDToString(List<int> aList, int timeStep, int targetTime,  bool hIgnored)
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
                int type = atomInputData.IdToType[ID - 1];
                string atomName = atomInputData.TypeToAtom[type];
                if (timeStep == targetTime)
                {
                    if (!IsListEquall(atomInputData.AtomList[timeStep][ID - 1], atomInputData.AtomList[timeStep + 1][ID - 1]))
                    {
                        if (atomName == "O" && atomInputData.AtomList[timeStep][ID - 1].Count == 1)
                        {
                            atomChain += "=";
                        }
                        atomChain += ("[" + atomName + ":1]");
                        continue;
                    }
                }
                else
                {
                    if (!IsListEquall(atomInputData.AtomList[timeStep][ID - 1], atomInputData.AtomList[targetTime][ID - 1]))
                    {
                        if (atomName == "O" && atomInputData.AtomList[timeStep][ID - 1].Count == 1)
                        {
                            atomChain += "=";
                        }
                        atomChain += ("[" + atomName + ":1]");
                        continue;
                    }
                }
                if (atomName == "O" && atomInputData.AtomList[timeStep][ID - 1].Count == 1)
                {
                    atomChain += "=";
                }
                if (atomName == "H" && hIgnored)
                {
                    atomChain += "X";
                    continue;
                }
                atomChain += atomName;
            }
            return atomChain;
        }

        static public bool IsListEquall(List<int> a, List<int> b)
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
    }
}
