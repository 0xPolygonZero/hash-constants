#
# This file is based on the Poseidon implementation at
# https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/659de89cd207e19b92852458dce92adf83ad7cf7/code/poseidonperm_x3_64_24_optimized.sage
#

# These come from ChaCha initialised with 0.
ROUND_CONSTANTS = [
    0xb585f767417ee042, 0x7746a55f77c10331, 0xb2fb0d321d356f7a, 0x0f6760a486f1621f,
    0xe10d6666b36abcdf, 0x8cae14cb455cc50b, 0xd438539cf2cee334, 0xef781c7d4c1fd8b4,
    0xcdc4a23a0aca4b1f, 0x277fa208d07b52e3, 0xe17653a300493d38, 0xc54302f27c287dc1,
    0x8628782231d47d10, 0x59cd1a8a690b49f2, 0xc3b919ad9efec0b0, 0xa484c4c637641d97,
    0x308bbd23f191398b, 0x6e4a40c1bf713cf1, 0x9a2eedb7510414fb, 0xe360c6e111c2c63b,
    0xd5c771901d4d89aa, 0xc35eae076e7d6b2f, 0x849c2656d0a09cad, 0xc0572c8c5cf1df2b,
    0xe9fa634a883b8bf3, 0xf56f6d4900fb1fdd, 0xf7d713e872a72a1b, 0x8297132b6ba47612,
    0xad6805e12ee8af1c, 0xac51d9f6485c22b9, 0x502ad7dc3bd56bf8, 0x57a1550c3761c577,
    0x66bbd30e99d311da, 0x0da2abef5e948f87, 0xf0612750443f8e94, 0x28b8ec3afb937d8c,
    0x92a756e6be54ca18, 0x70e741ec304e925d, 0x019d5ee2b037c59f, 0x6f6f2ed7a30707d1,
    0x7cf416d01e8c169c, 0x61df517bb17617df, 0x85dc499b4c67dbaa, 0x4b959b48dad27b23,
    0xe8be3e5e0dd779a0, 0xf5c0bc1e525ed8e6, 0x40b12cbf263cf853, 0xa637093f13e2ea3c,
    0x3cc3f89232e3b0c8, 0x2e479dc16bfe86c0, 0x6f49de07d6d39469, 0x213ce7beecc232de,
    0x5b043134851fc00a, 0xa2de45784a861506, 0x7103aaf97bed8dd5, 0x5326fc0dbb88a147,
    0xa9ceb750364cb77a, 0x27f8ec88cc9e991f, 0xfceb4fda8c93fb83, 0xfac6ff13b45b260e,
    0x7131aa455813380b, 0x93510360d5d68119, 0xad535b24fb96e3db, 0x4627f5c6b7efc045,
    0x645cf794e4da78a9, 0x241c70ed1ac2877f, 0xacb8e076b009e825, 0x3737e9db6477bd9d,
    0xe7ea5e344cd688ed, 0x90dee4a009214640, 0xd1b1edf7c77e74af, 0x0b65481bab42158e,
    0x99ad1aab4b4fe3e7, 0x438a7c91f1a360cd, 0xb60de3bd159088bf, 0xc99cab6b47a3e3bb,
    0x69a5ed92d5677cef, 0x5e7b329c482a9396, 0x5fc0ac0829f893c9, 0x32db82924fb757ea,
    0x0ade699c5cf24145, 0x7cc5583b46d7b5bb, 0x85df9ed31bf8abcb, 0x6604df501ad4de64,
    0xeb84f60941611aec, 0xda60883523989bd4, 0x8f97fe40bf3470bf, 0xa93f485ce0ff2b32,
    0x6704e8eebc2afb4b, 0xcee3e9ac788ad755, 0x510d0e66062a270d, 0xf6323f48d74634a0,
    0x0b508cdf04990c90, 0xf241708a4ef7ddf9, 0x60e75c28bb368f82, 0xa6217d8c3f0f9989,
    0x7159cd30f5435b53, 0x839b4e8fe97ec79f, 0x0d3f3e5e885db625, 0x8f7d83be1daea54b,
    0x780f22441e8dbc04, 0xeb9158465aedacd3, 0xd19e120d826c1b6c, 0x016ee53a7f007110,
    0xcb5fd54ed22dd1ca, 0xacb84178c58de144, 0x9c22190c2c463227, 0x5d693c1bcc98406d,
    0xdcef0798235f321a, 0x3d639263f55e0b1e, 0xe273fd977edb8fda, 0x418f027049d10fe7,
    0x8c25fda3f253a284, 0x2cbaed4dc25a884e, 0x5f58e6aff78dc2af, 0x284650ac6fb9d206,
    0x635b337f1391c13c, 0x9f9a036f1ac6361f, 0xb93e260cff6747b4, 0xb0a7eae8c7272e33,
    0xd0762cbce7da0a9f, 0x34c6efb829c754d6, 0x40bf0ab6166855c1, 0xb6b570fccc46a242,
    0x5a27b90055549545, 0xb1a5b166048b306f, 0x8722e0ad24f1006d, 0x788ee3b3b315049a,
    0x14a726661e5b0351, 0x98b7672fe1c3f13e, 0xbb93ae77bdc3aa8f, 0x28fd3b04756fc222,
    0x30a46805a86d7109, 0x337dc00c7844a0e7, 0xd5eca245253c861b, 0x77626382990d8546,
    0xc1e434bf33c3ae7a, 0x0299351a54dbf35e, 0xb2d456e4fb620184, 0x3e9ed1fdc00265ea,
    0x2972a92bb672e8db, 0x20216dd789f333ec, 0xadffe8cf746494a1, 0x1c4dbb1c5889d420,
    0x15a16a8a8c9972f5, 0x388a128b98960e26, 0x2300e5d6ca3e5589, 0x2f63aa865c9ceb9f,
    0xf1c36ce8d894420f, 0x271811252953f84a, 0xe5840293d5466a8e, 0x4d9bbc3e24e5f20e,
    0xea35bc29cfa2794b, 0x18e21b4bf59e2d28, 0x1e3b9fc632ef6adb, 0x25d643627a05e678,
    0x5a3f1bb1ecb63263, 0xdb7f0238ca031e31, 0xb462065960bfc4c4, 0x49c24ae463c280f4,
    0xd793862c6f7b901a, 0xaadd1106bdce475e, 0xc43b6e0eed8ad58f, 0xe29024c1f2060cb7,
    0x5e50c2755efbe17a, 0x10383f20ac183625, 0x38e8ee9d8a8a435d, 0xdd511837bcc52452,
    0x7750059861a7da6a, 0x86ab99b518d1dbef, 0xb1204f608ccfe33b, 0xef61ac84d8dfca49,
    0x1bbcd90f1f4eff36, 0x0cd1dabd9be9850a, 0x11a3ae5bf354bb11, 0xf755bfef11bb5516,
    0xa3b832506e2f3adb, 0x516306f4b617e6ba, 0xddb4ac4a2aeead3a, 0x64bb6dec62af4430,
    0xf9cc95c29895a152, 0x08d37f75632771b9, 0xeec49b619cee6b56, 0xf143933b56b3711a,
    0xe4c5dd82b9f6570c, 0xe7ad775756eefdc4, 0x92c2318bc834ef78, 0x739c25f93007aa0a,
    0x5636caca1725f788, 0xdd8f909af47cd0b6, 0xc6401fe16bc24d4e, 0x8ad97b342e6b3a3c,
    0x0c49366bb7be8ce2, 0x0784d3d2f4b39fb5, 0x530fb67ec5d77a58, 0x41049229b8221f3b,
    0x139542347cb606a3, 0x9cb0bd5ee62e6438, 0x02e3f615c4d3054a, 0x985d4f4adefb64a0,
    0x775b9feb32053cde, 0x304265a64d6c1ba6, 0x593664c3be7acd42, 0x4f0a2e5fd2bd6718,
    0xdd611f10619bf1da, 0xd8185f9b3e74f9a4, 0xef87139d126ec3b3, 0x3ba71336dd67f99b,
    0x7d3a455d8d808091, 0x660d32e15cbdecc7, 0x297a863f5af2b9ff, 0x90e0a736e6b434df,
    0x549f80ce7a12182e, 0x0f73b29235fb5b84, 0x16bf1f74056e3a01, 0x6d1f5a593019a39f,
    0x02ff876fa73f6305, 0xc5cb72a2fb9a5bd7, 0x8470f39d674dfaa3, 0x25abb3f1e41aea30,
    0x23eb8cc9c32951c7, 0xd687ba56242ac4ea, 0xda8d9e915d2de6b7, 0xe3cbdc7d938d8f1e,
    0xb9a8c9b4001efad6, 0xc0d28a5c64f2285c, 0x45d7ac9b878575b8, 0xeeb76e39d8da283e,
    0x3d06c8bd2fc7daac, 0x9c9c9820c13589f5, 0x65700b51db40bae3, 0x911f451579044242,
    0x7ae6849ff1fee8cc, 0x3bb340ebba896ae5, 0xb46e9d8bb71f0b4b, 0x8dcf22f9e1bde2a3,
    0x77bdaeda8cc55427, 0xf19e400ababa0e12, 0xc368a34939eb5c7f, 0x9ef1cd612c03bc5e,
    0xe89cd8553b94bbd8, 0x5cd377dcb4550713, 0xa7b0fb78cd4c5665, 0x7684403ef76c7128,
    0x5fa3f06f79c4f483, 0x8df57ac159dbade6, 0x2db01efa321b2625, 0x54846de4cfd58cb6,
    0xba674538aa20f5cd, 0x541d4963699f9777, 0xe9096784dadaa548, 0xdfe8992458bf85ff,
    0xece5a71e74a35593, 0x5ff98fd5ff1d14fd, 0x83e89419524c06e1, 0x5922040b6ef03286,
    0xf97d750eab002858, 0x5080d4c2dba7b3ec, 0xa7de115ba038b508, 0x6a9242acb5f37ec0,
    0xf7856ef865619ed0, 0x2265fc930dbd7a89, 0x17dfc8e5022c723b, 0x9001a64248f2d676,
    0x90004c13b0b8b50e, 0xb932b7cfc63485b0, 0xa0b1df81fd4c2bc5, 0x8ef1dd26b594c383,
    0x0541a4f9d20ba562, 0x9e611061be0a3c5b, 0xb3767e80e1e1624a, 0x0098d57820a88c6b,
    0x31d191cd71e01691, 0x410fefafbf90a57a, 0xbdf8f2433633aea8, 0x9e8cd55b9cc11c28,
    0xde122bec4acb869f, 0x4d001fd5b0b03314, 0xca66370067416209, 0x2f2339d6399888c6,
    0x6d1a7918f7c98a13, 0xdf9a493995f688f3, 0xebc2151f4ded22ca, 0x03cc2ba8a2bab82f,
    0xd341d03844ad9a9b, 0x387cb5d273ab3f58, 0xbba2515f74a7a221, 0x7248fe7737f37d9c,
    0x4d61e56a7437f6b9, 0x262e963c9e54bef8, 0x59e89b097477d296, 0x055d5b52b9e47452,
    0x82b27eb36e430708, 0xd30094caf3080f94, 0xcf5cb38227c2a3be, 0xfeed4db701262c7c,
    0x41703f5391dd0154, 0x5eeea9412666f57b, 0x4cd1f1b196abdbc4, 0x4a20358594b3662b,
    0x1478d361e4b47c26, 0x6f02dc0801d2c79f, 0x296a202eeb03c4b6, 0x2afd6799aec20c38,
    0x7acfd96f3050383d, 0x6798ba0c380dfdd3, 0x34c6f57b3de02c88, 0x5736e1baf82eb8a0,
    0x20057d2a0e58b8de, 0x3dea5bd5eb6e1404, 0x16e50d89874a6a98, 0x29bff3eccbfba19a,
    0x475cd3207974793c, 0x18a42105cde34cfa, 0x023e7414b0618331, 0x151471081b52594b,
    0xe4a3dff23bdeb0f3, 0x01a8d1a588c232ef, 0x11b4c74ee221d621, 0xe587cc0dce129c8c,
    0x1ff7327025a65080, 0x594e29c44b8602b1, 0xf6f31db1f5a56fd3, 0xc02ac5e4c7258a5e,
    0xe70201e9c5dc598f, 0x6f90ff3b9b3560b2, 0x42747a7262faf016, 0xd1f507e496927d26,
    0x1c86d265fdd24cd9, 0x3996ce73f6b5266e, 0x8e7fba02d68a061e, 0xba0dec71548b7546,
    0x9e9cbd785b8d8f40, 0xdae86459f6b3828c, 0xdebe08541314f71d, 0xa49229d29501358f,
    0x7be5ba0010c4df7c, 0xa3c95eaf09ecc39c, 0x0230bca8f5d457cd, 0x4135c2bedc68cdf9,
    0x166fc0cc4d5b20cc, 0x3762b59aa3236e6e, 0xe8928a4ceed163d2, 0x2a440b51b71223d9,
    0x80cefd2bb5f48e46, 0xbb9879c738328b71, 0x6e7c8f1ab47cced0, 0x164bb2de257ffc0a,
    0xf3c12fe5b800ea30, 0x40b9e92309e8c7e1, 0x551f5b0fe3b8d017, 0x25032aa7d4fc7aba,
    0xaaed340795de0a0a, 0x8ffd96bc38c8ba0f, 0x70fc91eb8aa58833, 0x7f795e2a97566d73,
    0x4543d9df72c4831d, 0xf172d73e69f20739, 0xdfd1c4ff1eb3d868, 0xbc8dfb62d26376f7,
    0x3f3b0b53ae4624f0, 0x8aff6a4012784bf9, 0xa788db8140374349, 0x7519463b9ca12e3f,
    0x916257fce69a385e, 0x1e2333a9f193f20a, 0x7e218f76de8e895d, 0x3a679f3e1277d39c,
    0x8291da1124f8da28, 0xa6915ec6c00589aa, 0xbd7f43cbef9842e6, 0x396046345ece5f81,
    0x9c0f589411f8f696, 0x97bf800069ca178d, 0xb70dadd0f6d9b373, 0xd3a8a413f5b0947a,
    0x710dd9eaee2506b7, 0x03a826cc0d051fb6, 0x78c126112bf8cb54, 0x5d13d03c9d6032e3,
    0x502a9452275bac54, 0xa399ed365a014c85, 0xd91d7776476759aa, 0x3388ea2010484f20,
    0x8bc8d8e0335c157b, 0x967830d0be15dfe6, 0x6b69b32c4e880930, 0xd730fbf87d6e70b5,
    0x85223b3b3b7bed07, 0xfa0d85221f3d72cf, 0xa6ded7c7f5d6d272, 0x9fbddba3c86edbad,
    0x0eac3c714b886b3d, 0x694d74d6b0a32710, 0x6b8ee6b99cd2ac11, 0x3501d69447d6f2b9,
    0x0290abe42804865e, 0x0331884444ffdab9, 0xf29c228bf96bc677, 0xf1d3df53724b6a6e,
    0xc23624ec8fd59daf, 0xb32c718470b115c2, 0x522741569915690a, 0xf4e786b5d6d87ecf,
    0x2c0d326a4e938885, 0xaff56278b67a71c2, 0x2ca8e42395fb3398, 0x8091f2239b333ebd,
    0xea4505d1a9901ff3, 0xdc48db966aac2a52, 0x48c7c93ff2349004, 0xa567c362fbb799da,
    0xd390081f9c257d4b, 0x9384c1070eb42745, 0xda4080d5e8af4bb3, 0x1929674becefc8e3,
    0xdce4ee6ff917b599, 0x613df420085e3c40, 0xe35ffbe87d486774, 0x19c8fe60f374b898,
    0x4ee06e1a5d0b2b59, 0xeb11f532da72e497, 0xa9f8b9b3b64cebcb, 0xd4e2607ab772baf5,
    0x082e4328a115b4bb, 0x571b0bc42e1d7204, 0xe5b2f2ab010ba6e8, 0xf95fca5d0b2e0dfd,
    0xa1515d97346e9a08, 0x25a98f6ef97ccefc, 0x0db6fa5ec2d790d0, 0x3df9019d524ffcf6,
    0x08ea1368d4262d01, 0xa09d948b08f421cd, 0x0ad5e966870973eb, 0x326ccef45586a421,
    0xc90f035d2f195b89, 0x17501e67f48218fc, 0xe18ee2416be34c37, 0xced93f36a645442c,
    0x754264af3522104f, 0x67eb06be66266452, 0xf37bfb01de53d2a6, 0x84dc8b4909daed76,
    0xfdab302f7750570f, 0x63e856e68d12f481, 0x957712a1717c53f8, 0x0eaf6bdd23ca146c,
    0x2c2f6e88585050a8, 0x1734b689a1096dbd, 0x435b37fe097ebbb4, 0xdb8ff81f76f13ea2,
    0xf26e6bf5ec48653b, 0xb3ef69ddf738923b, 0x0e5deaf37c2e7f7c, 0x537d810ff2ff5c7f,
    0x8b64439d81824d8f, 0x05cd0f2210bef923, 0x620236e39bcf6e15, 0x2f11cf1501e35f97,
    0x88752b7e6ff8173c, 0x2545e7e8883df01f, 0xdf7617a579ea6f69, 0x311b2fcdd633def4,
    0x179e1422f5d0c83f, 0x68d91273e79d1c5b, 0x3dd74da013b28e0e, 0x6d996a910a4a189d,
    0x5ebf824774f28933, 0x2b89e07840b58f52, 0x96397a640fa9b9c4, 0x90524d7159139759,
    0x258aa8b9125dabdb, 0x08f6bcce0819ce31, 0xe9469f8e3423540c, 0xde5164d3c42081e7,
    0x7d5c9ef70889cadd, 0xf9507747dd71c6c7, 0x121fd36677688512, 0x416e819aff900513,
    0x6196cce0fec4e80f, 0x3da35a4e648529ce, 0xa7f35230754bb327, 0x46f571b45c8489e1,
    0x8fd3c183a53fec29, 0xb040e8ca5ad907f6, 0xad672f5033fd6330, 0x771c671bf97c040f,
    0x460ab10909d9b3de, 0x216d1a78ce3a2bee, 0x6901012acef5593f, 0x7fd160fe46609bdb,
    0x0f37aece1b84d3aa, 0xe18d0ed7b26328cd, 0xfad907167e7a7c6c, 0x68a3629a69383872,
    0xeee91e2ffc681fea, 0xf6cd8d8714f06922, 0xd8cb08c311ec1ccd, 0xaf43fba5b4de25dd,
    0xba168eac7cabb1a3, 0xf45ea100a69a09b4, 0xdeffddfedfda58fc, 0xa64286d6cde2a535,
    0x96c7aeb47b963bdd, 0xac496870251ae7b3, 0x69e83510e83f6f00, 0x0af7b18fee7fa6bc,
    0xdc2fa9fbce94289c, 0x96050f80d386db64, 0x2ffc2568ebd20130, 0xeb131e5c04a334a0,
    0xbf1748518d9bc0b8, 0x93db9c179141eb5c, 0x34a9fb7bb2bce486, 0x0ca712203b21471f,
    0x9d4cf7a05a48207c, 0x7425b8e78f92e5ab, 0xb09a17910b15c6e6, 0x22c2c7adae292e9c,
    0x0c521e89ef7bfdfe, 0x6bf779019981aa32, 0xc691bb85d5b0bbb3, 0xeed836ab2c5012ec,
    0x38cd57fe93ee09a6, 0xe52a0af7fe66b04a, 0x06a3a9d191454450, 0xa8a147e1db6574ba,
    0xb832c1cb4c729ba0, 0x594acd3b9020d6da, 0xe0a8bea0820fd029, 0x7a73c734b55ef0aa,
    0x160d91cc51735394, 0x000728f94ce60a8d, 0x3091825206dce661, 0x9e235269a511546e,
    0xa7d9e03954f75d9b, 0x432f70e17905ede8, 0x472a2143e3799836, 0x0703c94d1dbc1d88,
    0x8d30fcff07a6eb4b, 0x9086dce770f1b9df, 0xaa6624bc8103a9cc, 0xcbd1b4ae92aac185,
    0x83332c17a793ee9b, 0x2bf7c5fc2d69b108, 0x1958860bcafd63b5, 0x9d49b02a68917c3c,
    0xdd956c968a83ba9c, 0xc809b28ff66f02fd, 0x64f7330c997c166d, 0x7e658282952e8b99,
    0x964addc37557be1f, 0x9294080ed3556c82, 0x934f3f9db6f863fe, 0xbf597ca5aa605d68,
    0x838e28334aa10ba5, 0x44311305ac3951bb, 0xa5d4622cdfcce3c1, 0xbf08aa3ec12c788c,
    0x0b09d598c90c7370, 0xa81b1af05019772d, 0x52ad553b02be9e9d, 0x956778e3b7b4eef6,
    0x351af23cc6e2c5ec, 0xf6ddb2c8918b7c20, 0x5a1b10b0dd4bc718, 0x6b573790e0190560,
    0x28f731229d0ea269, 0x6496d2b7f468d01a, 0xc37828a26e99fd51, 0x9f0d4d005eb173ce,
    0x97290a5cd121d82d, 0x733eb93c4523d1f5, 0xa485823c1dfaca30, 0x84a144f661077a8e,
    0xe0e2971d6e450e4d, 0xeb77ed68ba90c545, 0xfecb196cd5ccb5b0, 0x97348763528c945e,
    0xd0131ad33527cb5d, 0xd15b14a0e24318b1, 0xf809efb92933f86e, 0x93370f92e6ccfd30,
    0x77aacfbb4263fdc8, 0xd0209e3106220eee, 0xe3f45b2f7666b886, 0x1ee7aa194ac9e5fe,
    0x3a2a0ee62af22976, 0xc8096652ef57770c, 0x05e4ac783c617e58, 0x110882dca4469dd0,
    0xa734a6eb1f3af377, 0x95aa6b8cde4cc08b, 0x26170051bd7578cf, 0x61057afa25097f67,
    0x8f8a175ae3cd602c, 0x87a3ce5607bd038f, 0x4b0e54c0dc1d3718, 0xb1e7485936f01ba9,
]

class PoseidonData:
    def __init__(self, mds_matrix, sbox_exp, n_full_rounds, n_partial_rounds):
        assert mds_matrix.is_square()
        self.mds_matrix = mds_matrix
        self.state_width = self.mds_matrix.nrows()

        assert self.mds_matrix.base_ring().is_field(), 'matrix elements must be in a field'
        self.field = self.mds_matrix.base_ring()

        assert (self.field.cardinality() - 1) % sbox_exp != 0, 'sbox exponent must not divide p-1'
        self.sbox_exp = sbox_exp

        self.n_full_rounds = n_full_rounds
        self.n_partial_rounds = n_partial_rounds
        # This chunks the round_constants array into blocks of size state_width.
        self.round_constants = [vector(self.field, ROUND_CONSTANTS[index:index + self.state_width])
                                for index in range(0, len(ROUND_CONSTANTS), self.state_width)]
        self.round_constants_fast = calc_equivalent_constants(self)
        self.precomp_constants_fast = calc_equivalent_matrices(self)

#MDS_matrix = MDS_matrix.transpose() # QUICK FIX TO CHANGE MATRIX MUL ORDER (BOTH M AND M^T ARE SECURE HERE!)

def calc_equivalent_constants(hash_data):
    constants_temp = list(hash_data.round_constants)

    MDS_matrix_t = hash_data.mds_matrix.transpose()
    inv_MDS = MDS_matrix_t.inverse()

    # Start moving round constants up
    # Calculate c_i' = M^(-1) * c_(i+1)
    # Split c_i': Add c_i'[0] AFTER the S-box, add the rest to c_i
    # I.e.: Store c_i'[0] for each of the partial rounds, and make c_i = c_i + c_i' (where now c_i'[0] = 0)
    num_rounds = hash_data.n_full_rounds + hash_data.n_partial_rounds
    R_f = hash_data.n_full_rounds // 2
    for i in range(num_rounds - 2 - R_f, R_f - 1, -1):
        inv_cip1 = constants_temp[i+1] * inv_MDS
        constants_temp[i] = constants_temp[i] + vector([0] + inv_cip1[1:].list())
        constants_temp[i+1] = vector([inv_cip1[0]] + [0] * (hash_data.state_width - 1))

    return constants_temp

def calc_equivalent_matrices(hash_data):
    # Following idea: Split M into M' * M'', where M'' is "cheap" and
    # M' can move before the partial nonlinear layer
    #
    # The "previous" matrix layer is then M * M'. Due to the
    # construction of M', the M[0,0] and v values will be the same for
    # the new M' (and I also, obviously)
    #
    # Thus: Compute the matrices, store the w_hat and v values

    # HL: TODO: why is this transposed?
    MDS_matrix_t = hash_data.mds_matrix.transpose()

    w_hats = []
    vs = []

    M_mul = MDS_matrix_t
    M_i = matrix(hash_data.field, hash_data.state_width, hash_data.state_width) # == zero

    for i in range(hash_data.n_partial_rounds - 1, -1, -1):
        M_hat = M_mul.submatrix(1, 1)

        w = vector(M_mul[1:, 0])  # M_mul.submatrix(row=1, col=0, ncols=1)
        w_hat = M_hat.inverse() * w
        w_hats.append(w_hat)

        v = vector(M_mul[0, 1:]) # M_mul.submatrix(row=0, col=1, nrows=1)
        vs.append(v)

        # Generate new M_i, and multiplication M * M_i for "previous" round
        one = matrix([1])
        M_i = block_diagonal_matrix([one, M_hat])
        M_mul = MDS_matrix_t * M_i

    vs.reverse()
    w_hats.reverse()
    return [M_i, vs, w_hats, MDS_matrix_t[0, 0]]


# Computes s*A where s is the state row vector and A is the matrix
#
#    [ M_00  | v  ]
#    [ ------+--- ]
#    [ w_hat | Id ]
#
# M_00 is a scalar, v is 1x(t-1), w_hat is (t-1)x1 and Id is the
# (t-1)x(t-1) identity matrix.
def cheap_matrix_mul(state, v, w_hat, M_00):
    col = vector([M_00] + w_hat.list())  # leftmost column

    # rest = s0 * v + state[shift up by 1]
    # rest[j] = state dot column_{j+1} of A
    rest = state[0] * v + state[1:]

    # new s0 = [M_00 | w^] dot [state]
    return vector([col * state] + rest.list())

def full_rounds(state, round_ctr, hash_data):
    R_f = hash_data.n_full_rounds // 2
    for r in range(0, R_f):
        # Round constants, nonlinear layer, matrix multiplication
        state += hash_data.round_constants[round_ctr]
        for i in range(0, hash_data.state_width):
            state[i] = state[i]^hash_data.sbox_exp

        state = hash_data.mds_matrix * state
        round_ctr += 1
    return state, round_ctr

def partial_rounds_orig(state, round_ctr, hash_data):
    for r in range(0, hash_data.n_partial_rounds):
        # Round constants, nonlinear layer, matrix multiplication
        state += hash_data.round_constants[round_ctr]
        state[0] = state[0]^hash_data.sbox_exp
        state = hash_data.mds_matrix * state
        round_ctr += 1
    return state, round_ctr

def partial_rounds_fast(state, round_ctr, hash_data):
    M_i, vs, w_hats, M_00 = hash_data.precomp_constants_fast

    # Initial constants addition
    state += hash_data.round_constants_fast[round_ctr]

    # First full matrix multiplication
    state = state * M_i
    for r in range(0, hash_data.n_partial_rounds):
        # Round constants, nonlinear layer, matrix multiplication
        state[0] = state[0]^hash_data.sbox_exp

        # Moved constants addition
        if r < (hash_data.n_partial_rounds - 1):
            round_ctr += 1
            state[0] = state[0] + hash_data.round_constants_fast[round_ctr][0]
        # Optimized multiplication with cheap matrices
        state = cheap_matrix_mul(state, vs[r], w_hats[r], M_00)
    round_ctr += 1
    return state, round_ctr

def poseidon(init_state, hash_data):
    assert len(init_state) == hash_data.state_width, \
        f'expected initial state length of {hash_data.state_width} but it was {len(init_state)}'
    round_ctr = 0

    state = init_state
    state, round_ctr = full_rounds(state, round_ctr, hash_data)
    state, round_ctr = partial_rounds_fast(state, round_ctr, hash_data)
    state, round_ctr = full_rounds(state, round_ctr, hash_data)

    return state

def poseidon_original(init_state, hash_data):
    assert len(init_state) == hash_data.state_width, \
        f'expected initial state length of {hash_data.state_width} but it was {len(init_state)}'
    round_ctr = 0

    state = init_state
    state, round_ctr = full_rounds(state, round_ctr, hash_data)
    state, round_ctr = partial_rounds_orig(state, round_ctr, hash_data)
    state, round_ctr = full_rounds(state, round_ctr, hash_data)

    return state


def print_hex_vectlst(vs, indent=''):
    for v in vs:
        print(f'{indent}    [', end='')
        for i, a in enumerate(v):
            if i > 0 and (i % 4) == 0:
                print(f'\n{indent}     ', end='')
            print(f'0x{int(a):016x}, ', end='')
        print('],')


def test_consistency(hash_data):
    F = hash_data.field
    t = hash_data.state_width
    inputs = [
        vector(F, [0]*t),     # Test input [0, 0, ..., 0]
        vector(F, [F(-1)]*t), # Test input [p-1, ..., p-1]
        vector(F, range(t)),  # Test input [0, 1, ..., t-1]
        vector(F, [0xb69ed321abbeffbb, 0xfb496d8c39b64e42, 0x274f1cfbb925c789, 0x9e846d2b9a56b834,
                   0xc7f297c0d48bc3b6, 0xb859ab1e45850a0a, 0x3244fe3bcb1244cb, 0xb98e1cfa647575de])]

    for input_words in inputs:
        orig_output = poseidon_original(input_words, hash_data)
        fast_output = poseidon(input_words, hash_data)
        assert orig_output == fast_output

def print_fast_partial_consts(hash_data):
    round_consts = hash_data.round_constants_fast
    F = hash_data.field
    t = hash_data.state_width
    R_F = hash_data.n_full_rounds
    R_P = hash_data.n_partial_rounds

    indent = '    '
    R_f = R_F // 2
    print(f'\n{indent}const FAST_PARTIAL_FIRST_ROUND_CONSTANT: [u64; WIDTH]  = [', end='')
    cnt = 0
    for i in range(t):
        if (cnt % 4) == 0:
            print(f'\n{indent}    ', end='')
        val = int(round_consts[R_f][i])
        print(f'0x{val:016x}, ', end='')
        cnt += 1
    print(f'\n{indent}];\n')

    print(f'\n{indent}const FAST_PARTIAL_ROUND_CONSTANTS: [u64; N_PARTIAL_ROUNDS - 1]  = [', end='')
    cnt = 0
    for i, rc in enumerate(round_consts):
        # Skip the first four full rounds, and the first partial round;
        # stop at the end of the partial rounds.
        if i >= (R_f + 1) and i <= (R_F + R_P - (R_f + 1)):
            if (cnt % 4) == 0:
                print(f'\n{indent}    ', end='')
            val = int(rc[0])  # Only need first entry as the others are zero
            print(f'0x{val:016x}, ', end='')
            cnt += 1
    print(f'\n{indent}];\n')

    M_i, vs, w_hats, M_00 = hash_data.precomp_constants_fast

    print(f'{indent}const FAST_PARTIAL_ROUND_VS: [[u64; WIDTH - 1]; N_PARTIAL_ROUNDS] = [')
    print_hex_vectlst(vs, indent)
    print(f'\n{indent}];\n')

    print(f'{indent}const FAST_PARTIAL_ROUND_W_HATS: [[u64; WIDTH - 1]; N_PARTIAL_ROUNDS] = [')
    print_hex_vectlst(w_hats, indent)
    print(f'\n{indent}];\n')

    print(f'{indent}// NB: This is in COLUMN-major order to support cache-friendly pre-multiplication.')
    print(f'{indent}const FAST_PARTIAL_ROUND_INITIAL_MATRIX: [[u64; WIDTH - 1]; WIDTH - 1] = [')
    print_hex_vectlst(M_i.submatrix(1,1).columns(), indent)
    print(f'\n{indent}];\n')

if __name__ == "__main__":
    R_F = 8
    R_P = 22
    crandall_prime = 2^64 - 9 * 2^28 + 1
    crandall_field = GF(crandall_prime)
    sbox_exp = 7

    state_width = 8
    if state_width > 12:
        raise ValueError(state_width)

    crandall_binary_mds12 = matrix.circulant(
        vector(crandall_field, [1, 1, 2, 1, 8, 32, 2, 256, 4096, 8, 65536, 1024]))
    MDS_matrix = crandall_binary_mds12.submatrix(0, 0, state_width, state_width)
    hash_data = PoseidonData(MDS_matrix, sbox_exp, R_F, R_P)

    test_consistency(hash_data)
    print_fast_partial_consts(hash_data)
