using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace ReaxManager
{
    public class SortSmiles
    {
        private List<string> molNumTosmiles;
        private List<string> molNumToString;
        public SortSmiles(List<string> molNumToSmiles, List<string> molNumTostring)
        {
            this.molNumTosmiles = new List<string>(molNumToSmiles);
            this.molNumToString = new List<string>(molNumToString);
        }

        public void SortSpecies(ref int totalSpecies, Dictionary<string, int> stringToSpeciesNum, Dictionary<string, string> smilesPlusH)
        {
            List<string> smiles = new List<string>(molNumTosmiles);
            List<int> smilesCount = new();
            for (int i = 0; i < smiles.Count; i++)
            {
                smiles[i] = smiles[i].Replace("[", "");
                smiles[i] = smiles[i].Replace("]", "");
                smiles[i] = smiles[i].Replace(":", "");
                smiles[i] = smiles[i].Replace("1", "");
                smiles[i] = smiles[i].Replace("H", "");
                if (smiles[i] == "=O")
                {
                    continue;
                }
                string smiles_temp = smiles[i];
                smiles[i] = smiles[i].Replace("=", "");
                smiles[i] = smiles[i].Replace("(", "");
                smiles[i] = smiles[i].Replace(")", "");
                if (!smilesCount.Contains(smiles[i].Count()))
                {
                    smilesCount.Add(smiles[i].Count());
                    stringToSpeciesNum.Add(smiles_temp, 1);
                    smilesPlusH.Add(smiles_temp, molNumToString[i]);
                }
                else
                {
                    foreach (string key in stringToSpeciesNum.Keys)
                    {
                        if (smiles[i].Length == key.Replace("=", "").Replace("(", "").Replace(")", "").Length)
                        {
                            stringToSpeciesNum[key] += 1;
                        }
                    }
                }
            }
        }

        //それぞれの分子を特定して、数え上げていく
        public List<int> CountTargetMoleculeNumbers(List<string> targetMoleculeSmilesList)
        {
            List<string> smilesList = new List<string>(molNumTosmiles);
            List<int> molecularWeightList = SmilesToMolecularWeight(targetMoleculeSmilesList);
            List<int> moleculeNumbers = new List<int>();
            for (int i=0; i<molecularWeightList.Count; i++)
            {
                moleculeNumbers.Add(0);
            }
            
            foreach (string smiles in smilesList)
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
                smiles2 = smiles2.Replace("=", "");
                for(int i=0; i<molecularWeightList.Count; i++)
                {
                    if(smiles2.Length % molecularWeightList[i] == 0 && smiles2.Length >= molecularWeightList[i])
                    {
                        moleculeNumbers[i] = moleculeNumbers[i] + smiles2.Length/ moleculeNumbers[i];
                    }
                }                
            }
            return moleculeNumbers;
        }

        private List<int> SmilesToMolecularWeight(List<string> smilesList)
        {
            List<int> molecularWeightList = new List<int>();
            for(int i=0; i<smilesList.Count; i++)
            {
                smilesList[i] = smilesList[i].Replace("H", "");
                smilesList[i] = smilesList[i].Replace("[", "");
                smilesList[i] = smilesList[i].Replace("]", "");
                smilesList[i] = smilesList[i].Replace(":", "");
                smilesList[i] = smilesList[i].Replace("1", "");
                smilesList[i] = smilesList[i].Replace("(", "");
                smilesList[i] = smilesList[i].Replace(")", "");
                smilesList[i] = smilesList[i].Replace("=", "");
                molecularWeightList.Add(smilesList[i].Count());
            }
            return molecularWeightList;
        }
    }
}
