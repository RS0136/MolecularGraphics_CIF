import bpy
from bpy.props import StringProperty, BoolProperty
from Bio.PDB import MMCIFParser
import math

bl_info = {
    "name": "B-Mol",
    "author": "",
    "version": (0, 0),
    "blender": (2, 90, 0),
    "location": "3Dビューポート > Sidebar",
    "description": "タンパク質の立体構造を表示するアドオン",
    "warning": "",
    "support": "TESTING",
    "wiki_url": "",
    "tracker_url": "",
    "category": "3D View"}

#アミノ酸分類
nonpolar_AA = ["GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP"]
polaruncharged_AA = ["SER", "THR", "ASN", "GLN", "TYR", "CYS"]
polarcharged_AA = ["LYS", "ARG", "HIS", "ASP", "GLU"]
noncanonical_AA = ["MSE", "SEC", "PYL"]

#共有結合の距離
COVALENT_BOND_DISTANCE = 2

class SAMPLE27_OT_CreateObject(bpy.types.Operator):

    bl_idname = "object.sample27_create_object"
    bl_label = "分子グラフィクス"
    bl_description = "分子を表示します"
    bl_options = {"REGISTER", "UNDO"}

    #メニューを実行したときに呼ばれる関数
    def execute(self, context):

        #パスを取得
        path = bpy.context.scene.path_prop_str
        print(path)

        #各条件を取得
        main = bpy.context.scene.main_prop_bool
        side = bpy.context.scene.side_prop_bool
        ball = bpy.context.scene.ball_prop_bool
        stick = bpy.context.scene.stick_prop_bool
        nonpolar = bpy.context.scene.nonpolar_prop_bool
        polaruncharged = bpy.context.scene.polaruncharged_prop_bool
        polarcharged = bpy.context.scene.polarcharged_prop_bool
        noncanonical = bpy.context.scene.noncanonical_prop_bool

        #表示するアミノ酸
        global nonpolar_AA, polaruncharged_AA, polarcharged_AA, noncanonical_AA
        
        AA = []
        if nonpolar or polaruncharged or polarcharged or noncanonical:
            if nonpolar:
                AA += nonpolar_AA
            if polaruncharged:
                AA += polaruncharged_AA
            if polarcharged:
                AA += polarcharged_AA
            if noncanonical:
                AA += noncanonical_AA
        else:
            AA += nonpolar_AA + polaruncharged_AA + polarcharged_AA + noncanonical_AA
        
        #ファイルを開く
        cif_parser = MMCIFParser()
        structure = cif_parser.get_structure("", path)

        #原子座標として頂点，共有結合として辺を追加する
        global COVALENT_BOND_DISTANCE
        
        for model in structure.get_list():
            
            for chain in model.get_list():
                chain_atoms_coord = []
                chain_id = str(chain.get_id())
                
                for residue in chain.get_list():
                    if residue.get_resname() in AA:
                        
                        for atom in residue.get_list():
                            if main or side:
                                if main:
                                    if atom.get_name() in ["N", "CA", "C"]:
                                        chain_atoms_coord.append(atom.get_coord().tolist())
                                if side:
                                    if atom.get_name() not in ["N", "C"]:
                                        chain_atoms_coord.append(atom.get_coord().tolist())
                            else:
                                chain_atoms_coord.append(atom.get_coord().tolist())

                chain_bonds = []

                for a in range(len(chain_atoms_coord) - 1):
                    a_atom = chain_atoms_coord[a]
                    
                    for b in range(a + 1, len(chain_atoms_coord)):
                        b_atom = chain_atoms_coord[b]
                        
                        if ((a_atom[0] - b_atom[0]) ** 2 + (a_atom[1] - b_atom[1]) ** 2 + (a_atom[2] - b_atom[2]) ** 2) ** 0.5 < COVALENT_BOND_DISTANCE:
                            chain_bonds.append([a, b])

                if ball or stick:
                    if ball:
                        for atom in chain_atoms_coord:
                            bpy.ops.mesh.primitive_uv_sphere_add(radius = 1, location = atom)
                    if stick:
                        for atom in chain_atoms_coord:
                            bpy.ops.mesh.primitive_uv_sphere_add(radius = 0.5, location = atom)
                        for bond in chain_bonds:
                            a_atom = chain_atoms_coord[bond[0]]
                            b_atom = chain_atoms_coord[bond[1]]
                
                else:
                    msh = bpy.data.meshes.new(name = chain_id + "_mesh")
                    msh.from_pydata(chain_atoms_coord, chain_bonds, [])
                    obj = bpy.data.objects.new(name = chain_id, object_data = msh)
                    bpy.context.scene.collection.objects.link(obj)

        return {"FINISHED"}

#Sidebarのタブ[カスタムタブ]に，パネル[カスタムパネル]を追加
class SAMPLE27_PT_CustomPanel(bpy.types.Panel):

    bl_label = "分子グラフィクス" #パネルのヘッダに表示される文字列
    bl_space_type = "VIEW_3D" #パネルを登録するスペース
    bl_region_type = "UI" #パネルを登録するリージョン
    bl_category = "分子グラフィクス" #パネルを登録するタブ名
    bl_context = "objectmode" #パネルを表示するコンテキスト

    #ヘッダーのカスタマイズ
    def draw_header(self, context):
        layout = self.layout
        layout.label(text = "", icon = "PLUGIN")

    #メニューの描画処理
    def draw(self, context):
        layout = self.layout
        scene = context.scene

        #ファイルパスを入力するテキストボックスを追加
        layout.label(text = "PDBx/mmCIF File Path:")
        layout.prop(scene, "path_prop_str")

        #残基番号を入力するテキストボックスを追加
        layout.label(text = "Residue Number:")
        layout.prop(scene, "resnum_prop_str")

        #チェックボックスを追加
        layout.prop(scene, "main_prop_bool", text = "主鎖")
        layout.prop(scene, "side_prop_bool", text = "側鎖")
        layout.prop(scene, "ball_prop_bool", text = "ボール")
        layout.prop(scene, "stick_prop_bool", text = "スティック")
        layout.prop(scene, "nonpolar_prop_bool", text = "非極性")
        layout.prop(scene, "polaruncharged_prop_bool", text = "極性無電荷")
        layout.prop(scene, "polarcharged_prop_bool", text = "極性電荷")
        layout.prop(scene, "noncanonical_prop_bool", text = "Mse, Sec, Pyl")

        #ボタンを追加
        layout.operator(SAMPLE27_OT_CreateObject.bl_idname, text="表示")

#プロパティの初期化
def init_props():
    scene = bpy.types.Scene
    scene.path_prop_str = StringProperty(
        name = "パス",
        description = "PDBx/mmCIF file path",
        default = "")
    
    scene.resnum_prop_str = StringProperty(
        name = "残基番号",
        description = "Residue Number or Range",
        default = "")

    scene.main_prop_bool = BoolProperty(
        name = "主鎖",
        description = "Main on/off",
        default = False)

    scene.side_prop_bool = BoolProperty(
        name = "側鎖",
        description = "Side on/off",
        default = False)
    
    scene.ball_prop_bool = BoolProperty(
        name = "ボール",
        description = "Ball on/off",
        default = False)

    scene.stick_prop_bool = BoolProperty(
        name = "スティック",
        description = "Stick on/off",
        default = False)

    scene.nonpolar_prop_bool = BoolProperty(
        name = "非極性",
        description = "Non-Polar on/off",
        default = False)
    
    scene.polaruncharged_prop_bool = BoolProperty(
        name = "極性無電荷",
        description = "Polar Uncharged on/off",
        default = False)
    
    scene.polarcharged_prop_bool = BoolProperty(
        name = "極性電荷",
        description = "Polar Charged on/off",
        default = False)
    
    scene.noncanonical_prop_bool = BoolProperty(
        name = "Mse, Sec, Pyl",
        description = "Mse, Sec, Pyl on/off",
        default = False)

#プロパティを削除
def clear_props():
    scene = bpy.types.Scene
    del scene.path_prop_str, scene.resnum_prop_str
    del scene.main_prop_bool, scene.side_prop_bool
    del scene.ball_prop_bool, scene.stick_prop_bool
    del scene.nonpolar_prop_bool, scene.polaruncharged_prop_bool, scene.polarcharged_prop_bool, scene.noncanonical_prop_bool

classes = [
    SAMPLE27_PT_CustomPanel,
    SAMPLE27_OT_CreateObject]

def register():
    for c in classes:
        bpy.utils.register_class(c)
    init_props()
    print("アドオン「分子グラフィクス」が有効化されました。")

def unregister():
    clear_props()
    for c in classes:
        bpy.utils.unregister_class(c)
    print("アドオン「分子グラフィクス」が無効化されました。")

if __name__ == "__main__":
    register()
