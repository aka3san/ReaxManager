using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ReaxManager
{
    public class AtomInputData
    {
        private int totalTimeStep;
        public int TotalTimeStep => totalTimeStep;
        private int totalAtomNum;
        public int TotalAtomNum => totalAtomNum;
        private List<List<List<int>>> atomList;
        public List<List<List<int>>> AtomList => atomList;       
        private List<int> idToType;
        public List<int> IdToType => idToType;
        private Dictionary<int, string> typeToAtom;
        public Dictionary<int,string> TypeToAtom => typeToAtom;
        private List<string> targetMoleculeSmilesList;
        public List<string> TargetMoleculeSmilesList => targetMoleculeSmilesList;

        public AtomInputData(string filePath, int atomNum, int timeStep, Dictionary<int, string> typeToAtom, List<string> targetMoleculeSmilesList)
        {
            this.totalTimeStep = timeStep;
            this.atomList = new List<List<List<int>>>();
            this.targetMoleculeSmilesList = targetMoleculeSmilesList;
            for (int i = 0; i < timeStep; i++)
            {
                this.atomList.Add(new List<List<int>>());
                for (int j = 0; j < atomNum; j++)
                {
                    this.atomList[i].Add(new List<int>());
                }
            }
            this.idToType = new List<int>(atomNum);
            for (int i = 0; i < atomNum; i++)
            {
                this.idToType.Add(0);
            }
            this.typeToAtom = typeToAtom;       
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
                    this.atomList[timeStepCount][atomID - 1].Add(int.Parse(line[4 + j]));
                }
                this.idToType[atomID - 1] = int.Parse(line[2]);

            }
            IgnoreHydrogenBond();
            Ignore_C_();           
            this.totalAtomNum = this.atomList[0].Count;
        }

        public void IgnoreHydrogenBond()
        {
            for (int i = 0; i < totalTimeStep; i++)
            {
                for (int j = 0; j < atomList[i].Count; j++)
                {
                    if ((typeToAtom[idToType[j]] == "H" || typeToAtom[idToType[j]] == "F") && atomList[i][j].Count >= 2)
                    {
                        int mostChainedID = GetMostChainedAtom(i, j + 1);
                        for (int k = 0; k < atomList[i][j].Count; k++)
                        {
                            if (atomList[i][j][k] == mostChainedID)
                            {
                                continue;
                            }
                            atomList[i][atomList[i][j][k] - 1].Remove(j + 1);
                        }
                        atomList[i][j] = new List<int>() { mostChainedID };
                    }
                    atomList[i][j].RemoveAll(item => item == -1);
                }
            }
        }

        private void Ignore_C_()
        {
            for (int i = 0; i < totalTimeStep; i++)
            {
                for (int j = 0; j < atomList[i].Count; j++)
                {
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
        }

        public int GetMostChainedAtom(int timeStep, int atomID)
        {            
            int atom = this.atomList[timeStep][atomID - 1][0];
            foreach (int _atom in this.atomList[timeStep][atomID - 1])
            {
                if (this.atomList[timeStep][_atom - 1].Count > this.atomList[timeStep][atom - 1].Count)
                {
                    atom = _atom;
                }
            }
            return atom;            
        }

        public int RemoveMostChainedAtom(int timeStep, int atomID, List<List<int>> atomList_copy)
        {
            if (atomList_copy[atomID - 1] == null)
            {
                return -1;//結合している全ての原子を列挙済み
            }
            int atom = atomList_copy[atomID - 1][0];
            foreach (int _atom in atomList_copy[atomID - 1])
            {
                if (this.atomList[timeStep][_atom - 1].Count > this.atomList[timeStep][atom - 1].Count)
                {
                    atom = _atom;
                }
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
    }
    
}
