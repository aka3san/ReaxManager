using System;
using System.Collections.Generic;
namespace ReaxManager
{
    public class Main
    {
        public Main()
        {
            static void Main()
            {
                //タイプから原子の文字に変換するmap
                Dictionary<int, string> typeToAtom = new Dictionary<int, string>()
                {
                    {1,"O"},{ 2,"O"},{ 3,"H"},{ 4,"C"},{ 5,"C"},{ 6,"O"},{ 7,"C"},{ 8,"O"},
                    { 9,"F"},{ 10,"H"},{ 11,"H"},{ 12,"_C_"},{ 13,"_C_"},{ 14,"_C_"},{ 15,"H"},
                    { 16,"_C_"},{ 17,"_C_"},{ 18,"_C_"},{ 19, "_C_"}, { 20, "H"}, { 21, "_C_"}
                };

                // 対象となる分子(Smiles記法)の一覧  (これらの分子の分子数の経時変化が出力されます。分子数を解析したい分子を入力します。)
                List<string> targetMoleculeSmilesList = new List<string>()
                {
                    /*D-4Oh*/ "C(F)(F)(OC(F)(F)C(F)(F)C(F)(F)OC(F)(F)C(F)(F)C(F)(F)OC(F)(F)COCC(=O)CO)C(F)(F)C(F)(F)OC(F)(F)C(F)(F)C(F)(F)OC(F)(F)C(F)(F)C(F)(F)OC(F)(F)C(F)(F)C(F)(F)OC(F)(F)C(F)(F)C(F)(F)OC(F)(F)C(F)(F)C(F)(F)OC(F)(F)C(F)(F)C(F)(F)OC(F)(F)COCC(O)CO",
                    /*水*/ "O",
                    /*酸素*/ "O=O"
                };

                //この変数にインプット情報が格納されます。引数は左から順に、FilePath、総原子数(基盤も含める), 総タイムステップ数, typeから原子に変換するmap, 対象分子のリスト
                AtomInputData atomInputData = new AtomInputData("bondsH2OO2_short.reaxc", 48268, 252, typeToAtom, targetMoleculeSmilesList);

                //水素結合を無視する場合
                atomInputData.IgnoreHydrogenBond();

                //必要な情報を格納
                ReactionManager reactionManager = new ReactionManager(atomInputData);

                //テキストファイルに出力
                reactionManager.GetTextFileData();
            }
        }
    }
}
